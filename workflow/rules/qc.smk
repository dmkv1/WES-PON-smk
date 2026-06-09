rule fastqc:
    input:
        f"{config['outdir']}/bam/{{sample}}/{{sample}}.bam",
    output:
        html=f"{config['outdir']}/qc/fastqc/{{sample}}_fastqc.html",
    log:
        "logs/fastqc/fastqc_{sample}.log",
    conda:
        "../envs/qc.yaml"
    threads: 2
    params:
        out_dir=lambda wc, output: os.path.dirname(output.html),
    shell:
        "fastqc {input} -o {params.out_dir} -t {threads} > {log} 2>&1"


rule multiqc:
    input:
        fastp_html=expand(
            f"{config['outdir']}/qc/fastp/{{sample}}_fastp.html",
            sample=samples.index,
        ),
        fastqc_html=expand(
            f"{config['outdir']}/qc/fastqc/{{sample}}_fastqc.html",
            sample=samples.index,
        ),
        dupl_metrics=expand(
            f"{config['outdir']}/qc/metrics/{{sample}}.dupl_metrics.txt",
            sample=samples.index,
        ),
        mosdepth=expand(
            f"{config['outdir']}/qc/metrics/{{sample}}.mosdepth.summary.txt",
            sample=samples.index,
        ),
    output:
        f"{config['outdir']}/qc/multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/qc.yaml"
    params:
        outdir=config["outdir"],
    shell:
        "multiqc {params.outdir}/ logs/ -o {params.outdir}/qc --force > {log} 2>&1"
