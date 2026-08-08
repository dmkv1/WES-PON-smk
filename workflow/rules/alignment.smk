# FASTQ -> analysis-ready BAM, per GATK Best Practices ordering.
#
# This chain must stay command-for-command identical to WES-snakemake's
# workflow/rules/bam_mapping_gatk.smk. The normals built here are what every
# tumor in that pipeline is normalised against -- through the CNVkit reference,
# the Mutect2 PON and the PureCN normal database, all of which are functions of
# coverage, which duplicate marking moves. Two chains that differ produce a
# panel that quietly does not describe the samples it is applied to.
#
# The one deliberate omission is the caller's xengsort host-read filter, which
# runs for PDX samples only. A panel of normals contains no xenografts.


rule fastp_trim:
    input:
        fq1=lambda wc: get_fastq1(wc),
        fq2=lambda wc: get_fastq2(wc),
    output:
        fq1=temp("work/fastq/{sample}/{sample}.{unit}_R1.fq.gz"),
        fq2=temp("work/fastq/{sample}/{sample}.{unit}_R2.fq.gz"),
        html=f"{config['outdir']}/qc/fastp/{{sample}}.{{unit}}_fastp.html",
        json=f"{config['outdir']}/qc/fastp/{{sample}}.{{unit}}_fastp.json",
    log:
        "logs/fastp/fastp_{sample}.{unit}.log",
    conda:
        "../envs/fastp.yaml"
    threads: 4
    params:
        detect_adapter=(
            "--detect_adapter_for_pe"
            if config["params"]["fastp"]["detect_adapter_for_pe"]
            else ""
        ),
    shell:
        "fastp -i {input.fq1} -I {input.fq2} "
        "-o {output.fq1} -O {output.fq2} "
        "-h {output.html} -j {output.json} "
        "{params.detect_adapter} -w {threads} > {log} 2>&1"


# Per-unit alignment. Each unit is mapped with its own read group, resolved and
# validated at load time (workflow/scripts/units.py), so multi-lane samples
# carry distinct RGs and BQSR can model each unit's error profile separately. A
# single-FASTQ sample is just the one-unit case. Units are gathered by
# mark_duplicates below.
#
# The output is deliberately NOT coordinate-sorted: bwa emits query-grouped
# reads, samtools fixmate requires that grouping, and MarkDuplicates wants it
# too (see below). The single coordinate sort happens once, after duplicates
# are marked.
#
# -K fixes bwa's input chunk size so alignment is independent of thread count,
# which is what makes two runs of the same FASTQ comparable -- including runs of
# the two pipelines on hosts sized differently.
rule bwa_map:
    input:
        refg=config["refs"]["genome_human"],
        fq1="work/fastq/{sample}/{sample}.{unit}_R1.fq.gz",
        fq2="work/fastq/{sample}/{sample}.{unit}_R2.fq.gz",
    output:
        temp("work/bam/units/{sample}.{unit}.qgrp.bam"),
    params:
        rg=get_read_group,
    log:
        "logs/bwamem/bwamem_{sample}.{unit}.log",
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    resources:
        # The hg38 BWA index is ~5.5 GB resident and shared across threads; the
        # rest is per-thread alignment buffers plus the piped samtools fixmate.
        mem_mb=10240,
    shell:
        "(bwa mem -Y -K 100000000 -t {threads} -R '{params.rg}' "
        "{input.refg} {input.fq1} {input.fq2} "
        "| samtools fixmate -m -O bam,level=1 - {output}) 2> {log}"


# Gather all units of a sample and mark duplicates in one pass. Units sharing a
# library (LB) have their PCR duplicates detected across units; units from
# different libraries are treated separately, which is what preserves the
# independent evidence of two preps of one sample.
#
# ASSUME_SORT_ORDER queryname consumes bwa's query-grouped output directly, so
# secondary and supplementary reads are marked correctly. This is the ordering
# GATK Best Practices uses (MarkDuplicates before the coordinate sort), and it
# is why sort_bam below exists as its own rule. CREATE_INDEX is gone with it:
# a query-grouped BAM has no coordinate index to build.
rule mark_duplicates:
    input:
        get_unit_bams,
    output:
        bam=temp("work/bam/{sample}.md.qgrp.bam"),
        metrics=f"{config['outdir']}/qc/metrics/{{sample}}.dupl_metrics.txt",
    params:
        tmp_dir="tmp",
        inputs=lambda wildcards, input: " ".join(f"-I {b}" for b in input),
    log:
        "logs/MarkDuplicates/MarkDuplicates_{sample}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
        mem_mb=config["resources"]["mem_mb"],
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            MarkDuplicates \
            {params.inputs} \
            -O {output.bam} \
            -M {output.metrics} \
            --ASSUME_SORT_ORDER queryname \
            --TMP_DIR {params.tmp_dir} \
            >{log} 2>&1
        """


# samtools sort -m is a per-thread buffer, so the real footprint is
# sort_mem x threads. Without this the scheduler sees the 4000 MB default and
# can start a full-size GATK job beside the sort, which overcommits the host.
def _sort_mem_mb():
    raw = str(config["resources"]["sort_mem"]).strip()
    unit = raw[-1].upper()
    if unit in ("K", "M", "G"):
        value = float(raw[:-1])
        mb = {"K": value / 1024, "M": value, "G": value * 1024}[unit]
    else:
        mb = float(raw) / (1024 * 1024)  # a bare number is bytes, as in samtools
    # 15% above the buffers themselves, for the merge pass and the BGZF output.
    return int(mb * config["resources"]["threads"] * 1.15)


# The one coordinate sort. --write-index with the ##idx## form names the index
# {sample}.md.bai rather than {sample}.md.bam.bai, which is the filename the
# BQSR rules already expect.
rule sort_bam:
    input:
        "work/bam/{sample}.md.qgrp.bam",
    output:
        bam=temp("work/bam/{sample}.md.bam"),
        bai=temp("work/bam/{sample}.md.bai"),
    params:
        sort_mem=config["resources"]["sort_mem"],
        tmp_dir="tmp",
    log:
        "logs/sort_bam/sort_bam_{sample}.log",
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    resources:
        mem_mb=_sort_mem_mb(),
    shell:
        """
        mkdir -p {params.tmp_dir}
        samtools sort -@ {threads} -m {params.sort_mem} --write-index \
            -T {params.tmp_dir}/{wildcards.sample}.sort \
            -o {output.bam}##idx##{output.bai} {input} 2> {log}
        """
