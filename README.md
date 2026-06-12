# WES-PON-smk

A Snakemake pipeline that builds **Panels of Normals (PON)** for tumor-only Whole
Exome Sequencing (WES) analysis. Tumor-only calling needs a PON to filter
recurrent germline and artifact signals. From a set of normal samples this
pipeline produces, per capture-probe type:

- a **Mutect2 somatic PON** (`pon.vcf.gz`) for SNV/indel filtering, and
- **CNVkit reference profiles** (`reference_{sex}.cnn`, one per sex) for CNV calling.

It also aligns/recalibrates each normal and emits a MultiQC report.

## Workflow

![Pipeline rulegraph](rulegraph.png)

Samples are grouped by `probes` (capture kit) for the Mutect2 PON, and by
`probes × sex` for the CNVkit reference.

## Requirements

- Linux with **Singularity/Apptainer** (GATK and CNVkit run from pinned containers)
  and **conda** (fastp, bwa, samtools, mosdepth, fastqc, multiqc run from conda envs).
- Reference genome and resource files (see below) on local disk.

GATK 4.6.2.0 and CNVkit 0.9.13 images are pulled automatically; tool versions for
conda steps are pinned in `workflow/envs/`.

## Setup

### 1. Create the runner environment

```bash
conda env create -f environment.yml   # creates env "snakemake" (snakemake>=8.29.3, pandas)
conda activate snakemake
```

### 2. Provide reference files

Edit `config.yaml → refs` to point at your copies:

- `genome_human` — BWA-indexed reference FASTA (the Broad hg38 bundle, with `.fai`/`.dict`).
- `known_sites` — dbSNP, known indels, Mills (with `.idx`/`.tbi` siblings; GATK finds them).
- `refflat` — UCSC `refFlat.txt` (used by CNVkit autobin).

### 3. Declare your capture kits

Under `config.yaml → probe_configs`, each kit needs:

- `regions_bedfile` — the kit's **Regions** BED (GATK `--intervals`).
- `coverage_bedfile` — the kit's **Covered** BED, header-stripped (CNVkit target derivation).
- `library_prep` — label used in BAM read-group `LB`.

The key name (e.g. `SureSelectV6UTR`) is what you put in the `probes` column of the
sample sheet and becomes the `{probes}` wildcard in output paths.

### 4. Build the sample sheet

`samples.csv` with columns:

| column    | meaning |
|-----------|---------|
| `ID`      | unique sample id (used in all output filenames) |
| `sex`     | `m` or `f` (drives CNVkit `--sample-sex`; validated at startup) |
| `probes`  | a key from `probe_configs` |
| `R1`,`R2` | absolute paths to paired FASTQ files |

Read-group info (RGID/PU) is parsed from the first FASTQ header automatically.
Point to a different sheet via `config.yaml → samples_csv`.

### 5. Adjust run settings

- `profiles/default/config.yaml` — total `cores`, and the Singularity bind mount
  (`singularity-args`); add a `-B` for every host path your refs/FASTQs live under.
- `config.yaml → resources` — per-job `threads` and Java heap (`java_min_gb`/`java_max_gb`).
- Optional Telegram run notifications: set `telegram_bot_env` to a file exporting
  `TELEGRAM_BOT_TOKEN` and `TELEGRAM_CHAT_ID`; leave empty to disable. `run_id` labels messages.

## Running

```bash
./launch.sh           # full run; logs to snakemake.log
./launch.sh -n        # dry run (any extra args pass through to snakemake)
```

`launch.sh` invokes snakemake with `--profile profiles/default` (conda + singularity
enabled). Run from the pipeline root so relative paths (`work/`, `tmp/`, `logs/`) resolve.

## Outputs (`results/`)

- `PON/mutect2/{probes}/pon.vcf.gz` — Mutect2 somatic PON per probe type.
- `PON/cnvkit/{probes}/reference_{sex}.cnn` — CNVkit reference per probe type × sex.
- `bam/{sample}/{sample}.bam` — analysis-ready recalibrated BAMs.
- `qc/multiqc_report.html` — aggregated QC (fastp, FastQC, mosdepth, dup metrics).

Intermediates go to `work/`, scratch to `tmp/`, per-rule logs to `logs/`. All are
git-ignored.
