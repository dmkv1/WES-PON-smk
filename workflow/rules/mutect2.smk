rule mutect2_single_sample:
    input:
        bam=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bam",
        bai=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda wc: config["probe_configs"][probe_dict[wc.sample]][
            "regions_bedfile"
        ],
    output:
        vcf=f"{config['outdir']}/vcf/{{sample}}.vcf",
        idx=f"{config['outdir']}/vcf/{{sample}}.vcf.idx",
    log:
        "logs/Mutect2_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    threads: config["resources"]["threads"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            Mutect2 \
            -R {input.refg} -I {input.bam} \
            --intervals {input.regions} \
            --max-mnp-distance 0 \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """
