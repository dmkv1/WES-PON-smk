rule mosdepth:
    input:
        bam=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bam",
        bai=f"{config['outdir']}/{{sample}}/bam/{{sample}}.bai",
        regions_bed=lambda wc: config["probe_configs"][probe_dict[wc.sample]][
            "regions_bedfile"
        ],
    output:
        summary=f"{config['outdir']}/metrics/{{sample}}.mosdepth.summary.txt",
        region_dist=f"{config['outdir']}/metrics/{{sample}}.mosdepth.region.dist.txt",
        thresholds=f"{config['outdir']}/metrics/{{sample}}.thresholds.bed.gz",
    log:
        "logs/mosdepth_{sample}.log",
    conda:
        "../envs/qc.yaml"
    params:
        prefix=lambda wc, output: output.summary.replace(".mosdepth.summary.txt", ""),
    shell:
        "mosdepth --by {input.regions_bed} --thresholds 10,20,30,50 "
        "{params.prefix} {input.bam} > {log} 2>&1"
