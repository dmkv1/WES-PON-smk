rule purecn_coverage:
    # Reformats CNVkit's per-normal coverage into PureCN's legacy plain-text
    # format (Target/total_coverage/on_target columns). PureCN's own CNVkit
    # reader (.readCoverageCnn) never populates the `counts` field that
    # createNormalDatabase() requires internally, which crashes or silently
    # drops every interval — this format routes through .readCoverageGatk3
    # instead, which does populate it. Absolute scale of total_coverage
    # doesn't matter (PureCN only uses row-relative fractions), so
    # depth * width is a fine stand-in for a true read count.
    # Source is each normal's own .cnr (post cnvkit.py fix against the
    # sex-matched PON reference), not the raw target/antitarget .cnn pair:
    # tumors are fixed against that same reference downstream, and fix trims
    # bins with unusable signal. Feeding PureCN the untrimmed normal bin set
    # makes its identical() interval check reject every tumor comparison.
    input:
        cnr=f"{config['outdir']}/coverage/{{sample}}/{{sample}}.cnr",
    output:
        cov=f"{config['outdir']}/PON/purecn/coverage/{{sample}}.txt",
    run:
        cov = pd.read_csv(input.cnr, sep="\t")
        cov["on_target"] = cov["gene"] != "Antitarget"
        cov["Target"] = (
            cov["chromosome"]
            + ":"
            + (cov["start"] + 1).astype(str)
            + "-"
            + cov["end"].astype(str)
        )
        cov["total_coverage"] = cov["depth"] * (cov["end"] - cov["start"])
        cov[["Target", "total_coverage", "on_target"]].to_csv(
            output.cov, sep="\t", index=False
        )


rule purecn_coverage_list:
    # Split by sex, not just probe type: cnvkit_fix trims each normal's .cnr
    # to the bins usable in its sex-matched reference (reference_m.cnn vs
    # reference_f.cnn mask different bins, mostly chrX/chrY), so male and
    # female normals never share one interval set. PureCN's
    # createNormalDatabase() requires every coverage file in a NormalDB to
    # have identical intervals, so male and female normals must go into
    # separate NormalDBs.
    input:
        cov=lambda wc: expand(
            f"{config['outdir']}/PON/purecn/coverage/{{sample}}.txt",
            sample=samples[
                (samples["probes"] == wc.probes) & (samples["sex"] == wc.sex)
            ].index.tolist(),
        ),
    output:
        list=f"{config['outdir']}/PON/purecn/{{probes}}/coverage_files_{{sex}}.list",
    run:
        with open(output.list, "w") as f:
            f.write("\n".join(input.cov) + "\n")


rule haplotypecaller_gvcf:
    input:
        bam=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
        bai=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda wc: config["probe_configs"][probe_dict[wc.sample]][
            "covered_bedfile"
        ],
    output:
        gvcf=f"{config['outdir']}/gvcf/{{sample}}/{{sample}}.g.vcf.gz",
        tbi=f"{config['outdir']}/gvcf/{{sample}}/{{sample}}.g.vcf.gz.tbi",
    log:
        "logs/HaplotypeCaller/HaplotypeCaller_{sample}.log",
    container:
        config["containers"]["gatk"]
    threads: config["resources"]["threads"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            HaplotypeCaller \
            -R {input.refg} -I {input.bam} \
            -ERC GVCF \
            --intervals {input.regions} \
            -O {output.gvcf} \
            --native-pair-hmm-threads {threads} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """


rule genomicsdb_import_germline:
    # Separate GenomicsDB workspace from pon.smk's genomicsdb_import: that one
    # ingests Mutect2 tumor-only VCFs for the somatic PON, this one ingests
    # per-normal GVCFs for joint germline genotyping. Different objects, must
    # not share a workspace path.
    input:
        gvcfs=lambda wc: expand(
            f"{config['outdir']}/gvcf/{{sample}}/{{sample}}.g.vcf.gz",
            sample=samples[samples["probes"] == wc.probes].index.tolist(),
        ),
        tbis=lambda wc: expand(
            f"{config['outdir']}/gvcf/{{sample}}/{{sample}}.g.vcf.gz.tbi",
            sample=samples[samples["probes"] == wc.probes].index.tolist(),
        ),
        refg=config["refs"]["genome_human"],
        regions=lambda wc: config["probe_configs"][wc.probes]["covered_bedfile"],
    output:
        db=directory(f"work/genomicsdb_germline/{{probes}}"),
    log:
        "logs/GenomicsDBImport/GenomicsDBImportGermline_{probes}.log",
    container:
        config["containers"]["gatk"]
    threads: config["resources"]["threads"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        rm -rf {output.db}
        vcf_args=$(for v in {input.gvcfs}; do echo -V $v; done | tr '\\n' ' ')
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


rule genotype_gvcfs:
    input:
        db=f"work/genomicsdb_germline/{{probes}}",
        refg=config["refs"]["genome_human"],
    output:
        vcf=f"{config['outdir']}/PON/purecn/{{probes}}/normals_{{probes}}.joint.vcf.gz",
        tbi=f"{config['outdir']}/PON/purecn/{{probes}}/normals_{{probes}}.joint.vcf.gz.tbi",
    log:
        "logs/GenotypeGVCFs/GenotypeGVCFs_{probes}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            GenotypeGVCFs \
            -R {input.refg} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir} \
            >{log} 2>&1
        """


rule purecn_normaldb:
    # normal-panel stays pooled across both sexes per probe type: the
    # mapping-bias/beta-binomial step matches VCF samples by name against
    # the coverage list and tolerates a superset fine (confirmed — it
    # completed cleanly against the mixed-sex joint VCF before the
    # interval-mismatch crash below it). Only createNormalDatabase()'s
    # per-interval coverage matrix needs the sex split.
    input:
        coverage_list=f"{config['outdir']}/PON/purecn/{{probes}}/coverage_files_{{sex}}.list",
        joint_vcf=f"{config['outdir']}/PON/purecn/{{probes}}/normals_{{probes}}.joint.vcf.gz",
    output:
        normaldb=f"{config['outdir']}/PON/purecn/{{probes}}/normalDB_{{probes}}_{{sex}}_hg38.rds",
        mapping_bias=f"{config['outdir']}/PON/purecn/{{probes}}/mapping_bias_{{probes}}_{{sex}}_hg38.rds",
        hq_sites=f"{config['outdir']}/PON/purecn/{{probes}}/mapping_bias_hq_sites_{{probes}}_{{sex}}_hg38.bed",
        weights_png=f"{config['outdir']}/PON/purecn/{{probes}}/interval_weights_{{probes}}_{{sex}}_hg38.png",
        # NormalDB.R only writes this file when low.coverage.targets is
        # non-empty; touch it unconditionally afterward so Snakemake's
        # declared-output contract holds on kits where it isn't produced.
        low_coverage=f"{config['outdir']}/PON/purecn/{{probes}}/low_coverage_targets_{{probes}}_{{sex}}_hg38.bed",
    log:
        "logs/purecn_normaldb/purecn_normaldb_{probes}_{sex}.log",
    container:
        config["containers"]["purecn"]
    params:
        out_dir=lambda wc, output: os.path.dirname(output.normaldb),
        genome="hg38",
        assay=lambda wc: f"{wc.probes}_{wc.sex}",
    shell:
        """
        PURECN_SCRIPT=$(Rscript -e 'cat(system.file("extdata", "NormalDB.R", package="PureCN"))')
        Rscript $PURECN_SCRIPT \
            --out-dir {params.out_dir} \
            --genome {params.genome} \
            --assay {params.assay} \
            --coverage-files {input.coverage_list} \
            --normal-panel {input.joint_vcf} \
            --force \
            >{log} 2>&1
        touch {output.low_coverage}
        """
