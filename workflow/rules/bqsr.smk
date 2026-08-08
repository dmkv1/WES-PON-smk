rule base_recalibration:
    input:
        bam=f"work/bam/{{sample}}.md.bam",
        bai=f"work/bam/{{sample}}.md.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda wc: config["probe_configs"][probe_dict[wc.sample]][
            "covered_bedfile"
        ],
    output:
        recal_table=f"{config['outdir']}/qc/metrics/{{sample}}.recal_data.table",
    log:
        "logs/BaseRecalibrator/BaseRecalibrator_{sample}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
        mem_mb=config["resources"]["mem_mb"],
    params:
        tmp_dir="tmp",
        known_sites=lambda _: " ".join(
            f"--known-sites {s}" for s in config["refs"]["known_sites"]
        ),
        interval_padding=config["params"]["bqsr"]["interval_padding"],
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            BaseRecalibrator \
            -I {input.bam} -R {input.refg} \
            {params.known_sites} \
            --intervals {input.regions} \
            --interval-padding {params.interval_padding} \
            -O {output.recal_table} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """


rule apply_bqsr:
    input:
        bam=f"work/bam/{{sample}}.md.bam",
        recal_table=f"{config['outdir']}/qc/metrics/{{sample}}.recal_data.table",
        refg=config["refs"]["genome_human"],
    output:
        bam=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
        bai=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
    log:
        "logs/ApplyBQSR/ApplyBQSR_{sample}.log",
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
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            ApplyBQSR \
            -R {input.refg} -I {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.bam} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """
