SEX_MAP = {"m": "male", "f": "female"}


rule cnvkit_access:
    input:
        refg=config["refs"]["genome_human"],
    output:
        access=f"{config['outdir']}/cnvkit/access.bed",
    log:
        "logs/cnvkit_access.log",
    container:
        "docker://etal/cnvkit:0.9.11"
    threads: 1
    shell:
        "cnvkit.py access {input.refg} -o {output.access} > {log} 2>&1"


rule cnvkit_targets:
    input:
        coverage_bed=lambda wc: config["probe_configs"][wc.probes]["coverage_bedfile"],
        refflat=config["refs"]["refflat"],
        access=f"{config['outdir']}/cnvkit/access.bed",
    output:
        targets=f"{config['outdir']}/cnvkit/references/{{probes}}/targets.bed",
        antitargets=f"{config['outdir']}/cnvkit/references/{{probes}}/antitargets.bed",
    log:
        "logs/cnvkit_targets_{probes}.log",
    container:
        "docker://etal/cnvkit:0.9.11"
    threads: 1
    params:
        target_avg_size=config["params"]["cnvkit"]["target_avg_size"],
        antitarget_avg_size=config["params"]["cnvkit"]["antitarget_avg_size"],
    shell:
        """
        cnvkit.py target {input.coverage_bed} \
            --annotate {input.refflat} \
            --avg-size {params.target_avg_size} \
            -o {output.targets} >>{log} 2>&1
        cnvkit.py antitarget {output.targets} \
            -g {input.access} \
            --avg-size {params.antitarget_avg_size} \
            -o {output.antitargets} >>{log} 2>&1
        """


rule cnvkit_coverage:
    input:
        bam=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bam",
        bai=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bai",
        targets=lambda wc: (
            f"{config['outdir']}/cnvkit/references/"
            f"{probe_dict[wc.sample]}/targets.bed"
        ),
        antitargets=lambda wc: (
            f"{config['outdir']}/cnvkit/references/"
            f"{probe_dict[wc.sample]}/antitargets.bed"
        ),
    output:
        target_cov=f"{config['outdir']}/cnvkit/coverage/{{sample}}/{{sample}}.targetcoverage.cnn",
        antitarget_cov=f"{config['outdir']}/cnvkit/coverage/{{sample}}/{{sample}}.antitargetcoverage.cnn",
    log:
        "logs/cnvkit_coverage_{sample}.log",
    container:
        "docker://etal/cnvkit:0.9.11"
    threads: config["resources"]["threads"]
    shell:
        """
        cnvkit.py coverage {input.bam} {input.targets} \
            --count -p {threads} -o {output.target_cov} >>{log} 2>&1
        cnvkit.py coverage {input.bam} {input.antitargets} \
            --count -p {threads} -o {output.antitarget_cov} >>{log} 2>&1
        """


rule cnvkit_reference:
    input:
        target_covs=lambda wc: expand(
            f"{config['outdir']}/cnvkit/coverage/{{sample}}/{{sample}}.targetcoverage.cnn",
            sample=samples[
                (samples["probes"] == wc.probes) & (samples["sex"] == wc.sex)
            ].index.tolist(),
        ),
        antitarget_covs=lambda wc: expand(
            f"{config['outdir']}/cnvkit/coverage/{{sample}}/{{sample}}.antitargetcoverage.cnn",
            sample=samples[
                (samples["probes"] == wc.probes) & (samples["sex"] == wc.sex)
            ].index.tolist(),
        ),
        refg=config["refs"]["genome_human"],
    output:
        ref=f"{config['outdir']}/cnvkit/references/{{probes}}/reference_{{sex}}.cnn",
    log:
        "logs/cnvkit_reference_{probes}_{sex}.log",
    container:
        "docker://etal/cnvkit:0.9.11"
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
