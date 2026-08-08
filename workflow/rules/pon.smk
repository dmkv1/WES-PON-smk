rule genomicsdb_import:
    input:
        vcfs=lambda wc: expand(
            f"{config['outdir']}/vcf/{{sample}}/{{sample}}.vcf",
            sample=get_samples(wc.probes),
        ),
        idxs=lambda wc: expand(
            f"{config['outdir']}/vcf/{{sample}}/{{sample}}.vcf.idx",
            sample=get_samples(wc.probes),
        ),
        refg=config["refs"]["genome_human"],
        regions=lambda wc: config["probe_configs"][wc.probes]["covered_bedfile"],
    output:
        db=directory(f"work/genomicsdb/{{probes}}"),
    log:
        "logs/GenomicsDBImport/GenomicsDBImport_{probes}.log",
    container:
        config["containers"]["gatk"]
    threads: config["resources"]["threads"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
        mem_mb=config["resources"]["mem_mb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        rm -rf {output.db}
        vcf_args=$(for v in {input.vcfs}; do echo -V $v; done | tr '\\n' ' ')
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            GenomicsDBImport \
            -R {input.refg} -L {input.regions} \
            --merge-input-intervals \
            --batch-size 50 \
            --reader-threads {threads} \
            --genomicsdb-workspace-path {output.db} \
            $vcf_args \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """


rule create_somatic_pon:
    input:
        db=f"work/genomicsdb/{{probes}}",
        refg=config["refs"]["genome_human"],
    output:
        vcf=f"{config['outdir']}/PON/mutect2/{{probes}}/pon.vcf.gz",
        tbi=f"{config['outdir']}/PON/mutect2/{{probes}}/pon.vcf.gz.tbi",
    log:
        "logs/CreateSomaticPanelOfNormals/CreateSomaticPanelOfNormals_{probes}.log",
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
            CreateSomaticPanelOfNormals \
            -R {input.refg} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """
