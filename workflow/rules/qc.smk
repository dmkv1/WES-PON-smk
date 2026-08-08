# QC metric sources:
#   - FastQC reads the trimmed FASTQs (not the BQSR'd BAM, whose recalibrated
#     base qualities make FastQC's chemistry plots meaningless).
#   - CollectHsMetrics reports capture efficiency / on-target rate / fold
#     enrichment from the final recalibrated BAM, against the kit's bait+target
#     intervals (built once per probe kit by bed_to_interval_list).
#   - MultiQC depends on the raw metric files it actually parses (fastp .json,
#     fastqc .zip, Picard metrics, mosdepth summary), not the rendered HTML.

# Picard HsMetrics distinguishes baits (where probes hybridize -> the "Covered"
# BED) from targets (regions of interest -> the "Regions" BED).
_BED_FOR_KIND = {"bait": "covered_bedfile", "target": "target_regions_bedfile"}


rule bed_to_interval_list:
    input:
        bed=lambda wc: config["probe_configs"][wc.probes][_BED_FOR_KIND[wc.kind]],
        refg=config["refs"]["genome_human"],
    output:
        f"work/intervals/{{probes}}.{{kind}}.interval_list",
    log:
        "logs/BedToIntervalList/BedToIntervalList_{probes}_{kind}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
        mem_mb=config["resources"]["mem_mb"],
    wildcard_constraints:
        kind="bait|target",
    params:
        tmp_dir="tmp",
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            BedToIntervalList \
            -I {input.bed} -O {output} \
            -SD {input.refg} \
            -TMP_DIR {params.tmp_dir} \
            >{log} 2>&1
        """


rule fastqc:
    input:
        fq1="work/fastq/{sample}/{sample}.{unit}_R1.fq.gz",
        fq2="work/fastq/{sample}/{sample}.{unit}_R2.fq.gz",
    output:
        zip1=f"{config['outdir']}/qc/fastqc/{{sample}}.{{unit}}_R1_fastqc.zip",
        zip2=f"{config['outdir']}/qc/fastqc/{{sample}}.{{unit}}_R2_fastqc.zip",
    log:
        "logs/fastqc/fastqc_{sample}.{unit}.log",
    conda:
        "../envs/qc.yaml"
    threads: 2
    params:
        out_dir=lambda wc, output: os.path.dirname(output.zip1),
    shell:
        "fastqc {input.fq1} {input.fq2} -o {params.out_dir} -t {threads} > {log} 2>&1"


rule collect_hs_metrics:
    input:
        bam=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
        bai=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
        refg=config["refs"]["genome_human"],
        baits=lambda wc: f"work/intervals/{probe_dict[wc.sample]}.bait.interval_list",
        targets=lambda wc: f"work/intervals/{probe_dict[wc.sample]}.target.interval_list",
    output:
        metrics=f"{config['outdir']}/qc/metrics/{{sample}}.hs_metrics.txt",
    log:
        "logs/CollectHsMetrics/CollectHsMetrics_{sample}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
        mem_mb=config["resources"]["mem_mb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            CollectHsMetrics \
            -I {input.bam} -O {output.metrics} \
            -R {input.refg} \
            --BAIT_INTERVALS {input.baits} \
            --TARGET_INTERVALS {input.targets} \
            -TMP_DIR {params.tmp_dir} \
            >{log} 2>&1
        """


rule multiqc:
    input:
        # fastp and fastqc are per unit; everything computed on the finished BAM
        # is per sample. multiqc_config.yaml and the generated renames collapse
        # the two tiers back onto one row per sample where that is meaningful.
        fastp_json=[
            f"{config['outdir']}/qc/fastp/{sample}.{unit}_fastp.json"
            for sample in SAMPLES
            for unit in get_units(sample)
        ],
        fastqc_zip=[
            f"{config['outdir']}/qc/fastqc/{sample}.{unit}_R{read}_fastqc.zip"
            for sample in SAMPLES
            for unit in get_units(sample)
            for read in (1, 2)
        ],
        dupl_metrics=expand(
            f"{config['outdir']}/qc/metrics/{{sample}}.dupl_metrics.txt",
            sample=SAMPLES,
        ),
        mosdepth=expand(
            f"{config['outdir']}/qc/metrics/{{sample}}.mosdepth.summary.txt",
            sample=SAMPLES,
        ),
        hs_metrics=expand(
            f"{config['outdir']}/qc/metrics/{{sample}}.hs_metrics.txt",
            sample=SAMPLES,
        ),
        units_mqc=f"{config['outdir']}/qc/units_rg_mqc.tsv",
        row_type_mqc=f"{config['outdir']}/qc/row_type_mqc.tsv",
        config="multiqc_config.yaml",
        renames=f"{config['outdir']}/qc/multiqc_renames.yaml",
    output:
        f"{config['outdir']}/qc/multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/qc.yaml"
    params:
        outdir=config["outdir"],
    resources:
        # Parses every module's output for the whole cohort into one report.
        mem_mb=16384,
    shell:
        "multiqc -c {input.config} -c {input.renames} {params.outdir}/ logs/ "
        "-o {params.outdir}/qc --force > {log} 2>&1"
