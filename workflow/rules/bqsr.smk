rule base_recalibration:
    input:
        bam=f"{config['outdir']}/tmp/{{sample}}.md.bam",
        refg=config["refs"]["genome_human"],
    output:
        recal_table=f"{config['outdir']}/metrics/{{sample}}.recal_data.table",
    log:
        "logs/BaseRecalibrator_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
        known_sites=lambda _: " ".join(
            f"--known-sites {s}" for s in config["refs"]["known_sites"]
        ),
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            BaseRecalibrator \
            -I {input.bam} -R {input.refg} \
            {params.known_sites} \
            -O {output.recal_table} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """


rule apply_bqsr:
    input:
        bam=f"{config['outdir']}/tmp/{{sample}}.md.bam",
        recal_table=f"{config['outdir']}/metrics/{{sample}}.recal_data.table",
        refg=config["refs"]["genome_human"],
    output:
        bam=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bam",
        bai=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bai",
    log:
        "logs/ApplyBQSR_{sample}.log",
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
            ApplyBQSR \
            -R {input.refg} -I {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.bam} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """
