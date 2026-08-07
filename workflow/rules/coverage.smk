rule mosdepth:
    input:
        bam=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
        bai=f"{config['outdir']}/bam/{{sample}}/{{sample}}.bai",
        regions_bed=lambda wc: config["probe_configs"][probe_dict[wc.sample]][
            "covered_bedfile"
        ],
    output:
        summary=f"{config['outdir']}/qc/metrics/{{sample}}.mosdepth.summary.txt",
        region_dist=f"{config['outdir']}/qc/metrics/{{sample}}.mosdepth.region.dist.txt",
        thresholds=f"{config['outdir']}/qc/metrics/{{sample}}.thresholds.bed.gz",
    log:
        "logs/mosdepth/mosdepth_{sample}.log",
    conda:
        "../envs/qc.yaml"
    params:
        prefix=lambda wc, output: output.summary.replace(".mosdepth.summary.txt", ""),
    shell:
        "mosdepth --by {input.regions_bed} --thresholds 10,20,30,50 "
        "--no-per-base {params.prefix} {input.bam} > {log} 2>&1"
