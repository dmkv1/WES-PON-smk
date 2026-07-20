SEX_MAP = {"m": "male", "f": "female"}


rule cnvkit_access:
    input:
        refg=config["refs"]["genome_human"],
    output:
        access=f"{config['outdir']}/PON/cnvkit/access.bed",
    log:
        "logs/cnvkit_access/cnvkit_access.log",
    container:
        config["containers"]["cnvkit"]
    threads: 1
    shell:
        "cnvkit.py access {input.refg} -o {output.access} > {log} 2>&1"


rule cnvkit_strip_covered:
    # CNVkit chokes on the Agilent browser/track header lines, so strip them to a
    # plain BED here (GATK/mosdepth tolerate the header and use the file directly).
    input:
        covered=lambda wc: config["probe_configs"][wc.probes]["covered_bedfile"],
    output:
        temp(f"work/cnvkit/{{probes}}.covered.bed"),
    log:
        "logs/cnvkit_strip_covered/cnvkit_strip_covered_{probes}.log",
    shell:
        "grep -vE '^(browser|track|#)' {input.covered} > {output} 2> {log}"


rule cnvkit_autobin:
    input:
        covered_bed=f"work/cnvkit/{{probes}}.covered.bed",
        refflat=config["refs"]["refflat"],
        refg=config["refs"]["genome_human"],
        access=f"{config['outdir']}/PON/cnvkit/access.bed",
        bams=lambda wc: expand(
            f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
            sample=samples[samples["probes"] == wc.probes].index.tolist(),
        ),
        bais=lambda wc: expand(
            f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
            sample=samples[samples["probes"] == wc.probes].index.tolist(),
        ),
    output:
        targets=f"{config['outdir']}/PON/cnvkit/{{probes}}/targets.bed",
        antitargets=f"{config['outdir']}/PON/cnvkit/{{probes}}/antitargets.bed",
    log:
        "logs/cnvkit_autobin/cnvkit_autobin_{probes}.log",
    container:
        config["containers"]["cnvkit"]
    threads: 1
    params:
        method=config["params"]["cnvkit"]["autobin"]["method"],
        bp_per_bin=config["params"]["cnvkit"]["autobin"]["bp_per_bin"],
        target_min_size=config["params"]["cnvkit"]["autobin"]["target_min_size"],
        target_max_size=config["params"]["cnvkit"]["autobin"]["target_max_size"],
        antitarget_min_size=config["params"]["cnvkit"]["autobin"]["antitarget_min_size"],
        antitarget_max_size=config["params"]["cnvkit"]["autobin"]["antitarget_max_size"],
    shell:
        """
        cnvkit.py autobin {input.bams} \
            -m {params.method} \
            -t {input.covered_bed} \
            -g {input.access} \
            -f {input.refg} \
            --annotate {input.refflat} \
            -b {params.bp_per_bin} \
            --target-min-size {params.target_min_size} \
            --target-max-size {params.target_max_size} \
            --antitarget-min-size {params.antitarget_min_size} \
            --antitarget-max-size {params.antitarget_max_size} \
            --target-output-bed {output.targets} \
            --antitarget-output-bed {output.antitargets} \
            >{log} 2>&1
        """


rule cnvkit_coverage:
    input:
        bam=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
        bai=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
        targets=lambda wc: (
            f"{config['outdir']}/PON/cnvkit/" f"{probe_dict[wc.sample]}/targets.bed"
        ),
        antitargets=lambda wc: (
            f"{config['outdir']}/PON/cnvkit/"
            f"{probe_dict[wc.sample]}/antitargets.bed"
        ),
    output:
        target_cov=f"{config['outdir']}/coverage/{{sample}}/{{sample}}.targetcoverage.cnn",
        antitarget_cov=f"{config['outdir']}/coverage/{{sample}}/{{sample}}.antitargetcoverage.cnn",
    log:
        "logs/cnvkit_coverage/cnvkit_coverage_{sample}.log",
    container:
        config["containers"]["cnvkit"]
    threads: config["resources"]["threads"]
    shell:
        """
        cnvkit.py coverage {input.bam} {input.targets} \
            -p {threads} -o {output.target_cov} >>{log} 2>&1
        cnvkit.py coverage {input.bam} {input.antitargets} \
            -p {threads} -o {output.antitarget_cov} >>{log} 2>&1
        """


rule cnvkit_fix:
    # Normals need the same cnvkit.py fix pass tumors get downstream (fix
    # against the sex-matched PON reference), so their PureCN coverage is
    # built from the identical trimmed bin set instead of the raw,
    # untrimmed target+antitarget .cnn pair. Mismatched bin sets make
    # PureCN's identical() check reject every normal in the DB.
    input:
        target_cov=f"{config['outdir']}/coverage/{{sample}}/{{sample}}.targetcoverage.cnn",
        antitarget_cov=f"{config['outdir']}/coverage/{{sample}}/{{sample}}.antitargetcoverage.cnn",
        ref=lambda wc: (
            f"{config['outdir']}/PON/cnvkit/"
            f"{probe_dict[wc.sample]}/reference_{sex_dict[wc.sample]}.cnn"
        ),
    output:
        cnr=f"{config['outdir']}/coverage/{{sample}}/{{sample}}.cnr",
    log:
        "logs/cnvkit_fix/cnvkit_fix_{sample}.log",
    container:
        config["containers"]["cnvkit"]
    threads: 1
    shell:
        """
        cnvkit.py fix {input.target_cov} {input.antitarget_cov} {input.ref} \
            -o {output.cnr} >{log} 2>&1
        """


rule cnvkit_reference:
    input:
        target_covs=lambda wc: expand(
            f"{config['outdir']}/coverage/{{sample}}/{{sample}}.targetcoverage.cnn",
            sample=samples[
                (samples["probes"] == wc.probes) & (samples["sex"] == wc.sex)
            ].index.tolist(),
        ),
        antitarget_covs=lambda wc: expand(
            f"{config['outdir']}/coverage/{{sample}}/{{sample}}.antitargetcoverage.cnn",
            sample=samples[
                (samples["probes"] == wc.probes) & (samples["sex"] == wc.sex)
            ].index.tolist(),
        ),
        refg=config["refs"]["genome_human"],
    output:
        ref=f"{config['outdir']}/PON/cnvkit/{{probes}}/reference_{{sex}}.cnn",
    log:
        "logs/cnvkit_reference/cnvkit_reference_{probes}_{sex}.log",
    container:
        config["containers"]["cnvkit"]
    threads: config["resources"]["threads"]
    params:
        sex=lambda wc: SEX_MAP[wc.sex],
    shell:
        """
        cnvkit.py reference \
            {input.target_covs} {input.antitarget_covs} \
            --sample-sex {params.sex} \
            -f {input.refg} \
            -o {output.ref} >>{log} 2>&1
        """
