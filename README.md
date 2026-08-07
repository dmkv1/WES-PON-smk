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

`config.yaml` holds the run settings. The repository tracks `config.yaml.example` and
ignores `config.yaml`, so your host paths stay out of git. Copy the example first:

```bash
cp config.yaml.example config.yaml
```

Then repoint every `/path/to/...` placeholder. Under `refs`:

- `genome_human` — BWA-indexed reference FASTA (the Broad hg38 bundle, with `.fai`/`.dict`).
- `known_sites` — dbSNP, known indels, Mills (with `.idx`/`.tbi` siblings; GATK finds them).
- `refflat` — UCSC `refFlat.txt` (used by CNVkit autobin).

### 3. Declare your capture kits

Under `config.yaml → probe_configs`, each kit needs:

- `regions_bedfile` — the kit's **Regions** BED (GATK `--intervals`).
- `coverage_bedfile` — the kit's **Covered** BED, header-stripped (CNVkit target derivation).
- `library_prep` — label used in BAM read-group `LB`.

The key name (e.g. `SureSelectV6UTR`) becomes the `{probes}` wildcard in output paths,
and is what the calling pipeline's `panel_of_normals` paths point at. The sample sheet
does not carry it: each key declares a `capture_kit` token, and the sheet's
`capture_kit` column is matched against those tokens.

### 4. Build the sample sheet

`samplesheet.csv`, the same schema the calling pipeline uses:

| column           | meaning |
|------------------|---------|
| `ID`          | sample group identifier, usually a patient ID. No `_`, `.` or `/` |
| `sample`      | sample identifier, unique across the whole cohort (this pipeline has no run dimension to tell two same-named samples apart) |
| `gender`      | `m` or `f`. Drives the CNVkit `--sample-sex` and the per-sex PureCN normal database; validated at startup |
| `capture_kit` | a `capture_kit` token from `probe_configs` |
| `fq1`, `fq2`  | path to the R1/R2 FASTQ, or a glob matching several. A glob expands to one alignment unit per matched pair |

Optional `flowcell`, `lane`, `barcode` and `library` columns override the read groups
per row; otherwise they are derived from the FASTQ headers and filenames. The rules are
the caller's — `workflow/scripts/units.py` is vendored from it verbatim so that normals
and tumors resolve read groups identically.

The caller's schema adds `sample_type` and `tumor_fraction`. This pipeline takes no
position on either — every sample it is given is a normal, and there is no tumor
fraction to speak of — so both are accepted, carried through to
`results/metadata/samples.tsv` as provenance, and otherwise ignored. One generated
sheet therefore feeds both pipelines, and a sheet written for this one alone can leave
them out.

Point to a different sheet via `config.yaml → samplesheet`.

In this lab the sheet is generated from the catalogue:

```bash
python3 ../../scripts/ingest.py samplesheet --out samplesheet.csv --types CTRL,BL
```

Nothing in the pipeline reads that tool or its catalogue; any sheet with these columns
works, including a hand-written one.

### 5. Adjust run settings

Machine capacity lives in the run profile, which is ignored by git like `config.yaml`.
Copy its example too:

```bash
cp profiles/default/config.yaml.example profiles/default/config.yaml
```

Both examples are sized for a 16-thread, 64 GB machine. The two must be scaled
together, or the scheduler's ceiling stops matching what a job actually takes:

- `profiles/default/config.yaml` — total `cores` and the `resources.mem_mb` budget the
  scheduler packs against.
- `config.yaml → resources` — per-job `threads`, Java heap (`java_min_gb`/`java_max_gb`)
  and the `mem_mb` each GATK job declares.

The container bind mount needs no setting: the Snakefile derives it from `refs.path`,
so every reference and BED path must resolve under that root.
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
