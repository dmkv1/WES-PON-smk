rule mutect2_single_sample:
    input:
        bam=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
        bai=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda wc: config["probe_configs"][probe_dict[wc.sample]][
            "covered_bedfile"
        ],
    output:
        vcf=f"{config['outdir']}/vcf/{{sample}}/{{sample}}.vcf",
        idx=f"{config['outdir']}/vcf/{{sample}}/{{sample}}.vcf.idx",
        stats=f"{config['outdir']}/vcf/{{sample}}/{{sample}}.vcf.stats",
    log:
        "logs/Mutect2/Mutect2_{sample}.log",
    container:
        config["containers"]["gatk"]
    # pairHMM is the only threaded stage and it scales poorly past a few
    # threads. With one job per normal, throughput comes from running many
    # samples at once rather than widening each one.
    threads: 2
    resources:
        java_min_gb=config["resources"]["gatk"]["medium"]["java_min_gb"],
        java_max_gb=config["resources"]["gatk"]["medium"]["java_max_gb"],
        mem_mb=config["resources"]["gatk"]["medium"]["mem_mb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            Mutect2 \
            -R {input.refg} -I {input.bam} \
            --intervals {input.regions} \
            --max-mnp-distance 0 \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """
