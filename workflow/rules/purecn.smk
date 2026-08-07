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
        cov=f"{config['outdir']}/purecn/coverage/{{sample}}.txt",
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
            f"{config['outdir']}/purecn/coverage/{{sample}}.txt",
            sample=get_samples(wc.probes, wc.sex),
        ),
    output:
        list=f"{config['outdir']}/purecn/{{probes}}/coverage_files_{{sex}}.list",
    run:
        with open(output.list, "w") as f:
            f.write("\n".join(input.cov) + "\n")


CANONICAL_CONTIGS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


rule purecn_interval_check:
    # Interval-consistency gate for the CNV/purity arm. The PON's reference and
    # NormalDB must share one canonical interval set by construction, or PureCN's
    # identical() check rejects every tumor .cnr downstream. A prior verify pass
    # only checked sample names, not intervals — that gap let a 6-bin ALT-contig
    # divergence (chr22_KI270879v1_alt / GSTT1) reach production. This asserts:
    #   1. no non-canonical contig survives in targets/antitargets/reference, and
    #   2. every normal's coverage file feeding a NormalDB shares one interval set
    #      (the exact precondition createNormalDatabase() enforces internally).
    input:
        targets=f"{config['outdir']}/PON/cnvkit/{{probes}}/targets.bed",
        antitargets=f"{config['outdir']}/PON/cnvkit/{{probes}}/antitargets.bed",
        reference=f"{config['outdir']}/PON/cnvkit/{{probes}}/reference_{{sex}}.cnn",
        coverage_list=f"{config['outdir']}/purecn/{{probes}}/coverage_files_{{sex}}.list",
        normaldb=f"{config['outdir']}/PON/purecn/{{probes}}/normalDB_{{probes}}_{{sex}}_hg38.rds",
    output:
        ok=f"{config['outdir']}/purecn/{{probes}}/interval_check_{{sex}}.ok",
    run:
        def contigs(path, chrom_col, sep, header):
            # antitargets.bed is empty when off-target is disabled
            # (params.cnvkit.use_offtarget=false); pandas raises on a 0-byte
            # read, so treat it as the empty contig set.
            if os.path.getsize(path) == 0:
                return set()
            df = pd.read_csv(path, sep=sep, header=header)
            return set(df.iloc[:, chrom_col].astype(str))

        # 1. canonical-only across the CNVkit interval artifacts
        for name, path, col, sep, hdr in [
            ("targets", input.targets, 0, "\t", None),
            ("antitargets", input.antitargets, 0, "\t", None),
            ("reference", input.reference, 0, "\t", 0),
        ]:
            noncanon = contigs(path, col, sep, hdr) - CANONICAL_CONTIGS
            if noncanon:
                raise ValueError(
                    f"{name} ({path}) contains non-canonical contigs: "
                    f"{sorted(noncanon)}"
                )

        # 2. identical interval set across every normal in this NormalDB
        with open(input.coverage_list) as f:
            cov_files = [ln.strip() for ln in f if ln.strip()]
        ref_set = None
        for cov in cov_files:
            targets = tuple(pd.read_csv(cov, sep="\t")["Target"])
            if ref_set is None:
                ref_set, ref_file = targets, cov
            elif targets != ref_set:
                raise ValueError(
                    f"interval mismatch: {cov} ({len(targets)} bins) vs "
                    f"{ref_file} ({len(ref_set)} bins) — NormalDB would be "
                    f"rejected by createNormalDatabase()"
                )
        with open(output.ok, "w") as f:
            f.write(
                f"OK {wildcards.probes} {wildcards.sex}: "
                f"{len(cov_files)} normals, {len(ref_set)} canonical bins\n"
            )


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
            sample=get_samples(wc.probes),
        ),
        tbis=lambda wc: expand(
            f"{config['outdir']}/gvcf/{{sample}}/{{sample}}.g.vcf.gz.tbi",
            sample=get_samples(wc.probes),
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
        vcf=f"{config['outdir']}/purecn/{{probes}}/normals_{{probes}}.joint.vcf.gz",
        tbi=f"{config['outdir']}/purecn/{{probes}}/normals_{{probes}}.joint.vcf.gz.tbi",
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
        coverage_list=f"{config['outdir']}/purecn/{{probes}}/coverage_files_{{sex}}.list",
        joint_vcf=f"{config['outdir']}/purecn/{{probes}}/normals_{{probes}}.joint.vcf.gz",
    output:
        # Deliverables consumed by the calling pipeline stay under PON/.
        normaldb=f"{config['outdir']}/PON/purecn/{{probes}}/normalDB_{{probes}}_{{sex}}_hg38.rds",
        mapping_bias=f"{config['outdir']}/PON/purecn/{{probes}}/mapping_bias_{{probes}}_{{sex}}_hg38.rds",
        # NormalDB.R byproducts are internal build artifacts, kept out of PON/.
        hq_sites=f"{config['outdir']}/purecn/{{probes}}/mapping_bias_hq_sites_{{probes}}_{{sex}}_hg38.bed",
        weights_png=f"{config['outdir']}/purecn/{{probes}}/interval_weights_{{probes}}_{{sex}}_hg38.png",
        # NormalDB.R only writes this file when low.coverage.targets is
        # non-empty; touch it unconditionally afterward so Snakemake's
        # declared-output contract holds on kits where it isn't produced.
        low_coverage=f"{config['outdir']}/purecn/{{probes}}/low_coverage_targets_{{probes}}_{{sex}}_hg38.bed",
    log:
        "logs/purecn_normaldb/purecn_normaldb_{probes}_{sex}.log",
    container:
        config["containers"]["purecn"]
    params:
        # NormalDB.R writes every artifact into one --out-dir; point it at the
        # internal build dir, then move the two deliverables into PON/.
        out_dir=lambda wc, output: os.path.dirname(output.hq_sites),
        pon_dir=lambda wc, output: os.path.dirname(output.normaldb),
        genome="hg38",
        assay=lambda wc: f"{wc.probes}_{wc.sex}",
    shell:
        """
        mkdir -p {params.pon_dir}
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
        mv {params.out_dir}/$(basename {output.normaldb}) {output.normaldb}
        mv {params.out_dir}/$(basename {output.mapping_bias}) {output.mapping_bias}
        """
