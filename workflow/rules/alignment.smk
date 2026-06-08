rule fastp_trim:
    input:
        fq1=lambda wc: get_fastq1(wc),
        fq2=lambda wc: get_fastq2(wc),
    output:
        fq1=temp("work/fastq/{sample}/{sample}.trimmed.1.fq.gz"),
        fq2=temp("work/fastq/{sample}/{sample}.trimmed.2.fq.gz"),
        html=f"{config['outdir']}/qc/fastp/{{sample}}_fastp.html",
        json=f"{config['outdir']}/qc/fastp/{{sample}}_fastp.json",
    log:
        "logs/fastp_{sample}.log",
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


rule bwa_map:
    input:
        refg=config["refs"]["genome_human"],
        fq1="work/fastq/{sample}/{sample}.trimmed.1.fq.gz",
        fq2="work/fastq/{sample}/{sample}.trimmed.2.fq.gz",
    output:
        temp(f"{config['outdir']}/tmp/{{sample}}.raw.bam"),
    log:
        "logs/bwamem_{sample}.log",
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    shell:
        "(bwa mem -Y -t {threads} {input.refg} {input.fq1} {input.fq2} "
        "| samtools view -Sb - > {output}) 2> {log}"


rule add_read_groups:
    input:
        f"{config['outdir']}/tmp/{{sample}}.raw.bam",
    output:
        temp(f"{config['outdir']}/tmp/{{sample}}.rg.bam"),
    log:
        "logs/AddOrReplaceReadGroups_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        rg=get_read_group_params,
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            AddOrReplaceReadGroups \
            -I {input} -O {output} \
            -ID {params.rg[RGID]} -SM {params.rg[RGSM]} \
            -PU {params.rg[RGPU]} -LB {params.rg[RGLB]} \
            -PL {params.rg[RGPL]} -TMP_DIR {params.tmp_dir} \
            >{log} 2>&1
        """


rule fix_mate_info:
    input:
        f"{config['outdir']}/tmp/{{sample}}.rg.bam",
    output:
        temp(f"{config['outdir']}/tmp/{{sample}}.fixmate.bam"),
    log:
        "logs/FixMateInformation_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            FixMateInformation \
            -I {input} -O {output} \
            -SO coordinate -VALIDATION_STRINGENCY SILENT \
            -TMP_DIR {params.tmp_dir} \
            >{log} 2>&1
        """


rule mark_duplicates:
    input:
        f"{config['outdir']}/tmp/{{sample}}.fixmate.bam",
    output:
        bam=temp(f"{config['outdir']}/tmp/{{sample}}.md.bam"),
        metrics=f"{config['outdir']}/metrics/{{sample}}.dupl_metrics.txt",
    log:
        "logs/MarkDuplicates_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            MarkDuplicates \
            -I {input} -O {output.bam} -M {output.metrics} \
            --CREATE_INDEX false --TMP_DIR {params.tmp_dir} \
            >{log} 2>&1
        """
