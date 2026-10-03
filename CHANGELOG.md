# Changelog

This project uses [Semantic Versioning](https://semver.org), adapted for a pipeline.
The public contract is the samplesheet schema, the `config.yaml` keys and the output
paths under `results/PON/`, which the caller's `panel_of_normals` config points at.

Each release states whether it changes the PON artifacts for the same set of normals.
A version that changes them needs a rebuild before tumor-only calls are compared across
PON versions.

| Part | Reason |
|---|---|
| MAJOR | A samplesheet column, a `config.yaml` key or a deliverable path changes. A reference file becomes necessary. A deliverable is removed |
| MINOR | A tool, a rule or a deliverable is added, and an existing configuration still runs |
| PATCH | A bug fix, a documentation change, or a pinned version bump that does not change the artifacts |

## [1.2.0] - 2026-10-03

**With default settings this release does not change the PON artifacts.**
`config.yaml.example` sets two of the new options: `trim_front: 5` for V8+UTR, which
changes the V8+UTR reads and every V8+UTR artifact, and
`params.mutect2.interval_padding: 150`, which changes the Mutect2 PON of every kit.
Rebuild before use with WES-snakemake 2.4.0 or later configured the same way.

### Added
- `probe_configs.<kit>.trim_front` (default 0): fastp `--trim_front1/--trim_front2`,
  matching WES-snakemake 2.4.0. SureSelect XT HS2 libraries (V8+UTR) start each read
  with a 3-bp molecular barcode and 1-2 dark bases. Untrimmed, bwa soft-clips them in
  ~95% of reads.
- `params.mutect2.interval_padding` (default 0): `--interval-padding` for the normals'
  Mutect2 and for GenomicsDBImport, matching the caller's Mutect2 calling territory.

### Changed
- `profiles/default/config.yaml` is tracked and is a safe floor: 8 cores, 64 GB,
  `io_heavy: 2`. `profiles/default/config.yaml.example` is removed. Size a real run
  with a gitignored `profiles/<name>/` and `./launch.sh --workflow-profile <name>`.
  Move an existing untracked `profiles/default/config.yaml` to `profiles/<name>/`
  before pulling.
- `launch.sh` passes its arguments straight to snakemake and no longer adds
  `--profile profiles/default`, matching WES-snakemake 2.4.1.

## [1.1.0] - 2026-09-22

**Changes the PON artifacts. Rebuild before use.**

* `reference_m.cnn` is built with `--male-reference`. A normal male chrX sits at
  log2 0, which matches the caller's `cnvkit.py call --male-reference`. The male
  `.cnr` files and the male PureCN NormalDB change with it. `reference_f.cnn` is
  unchanged.
* `interval_check_{sex}.ok` also asserts that the median chrX log2 of every
  normal's `.cnr` lies within ±0.3 of 0.
* `units.py` is re-vendored from the caller. It accepts the optional
  `known_ploidy` column, which this pipeline ignores. Read-group resolution is
  unchanged.

## [1.0.0] - 2026-08-08

The first tagged release.

The pipeline aligns and recalibrates a set of normal samples and builds four
deliverables under `results/PON/`:

* the Mutect2 somatic panel of normals, per capture kit,
* the CNVkit bin definitions, per capture kit,
* the CNVkit reference profiles, per capture kit and sex,
* the PureCN normal database and mapping bias, per capture kit and sex.

It also writes a MultiQC report over the normals.

The preprocessing arm and the read-group resolution are vendored from
[WES-snakemake](https://github.com/dmkv1/WES-snakemake), so that a PON is built with the
same alignment and recalibration as the calls that consume it. `tests/test_vendored_units.py`
compares the vendored files against the caller's published ref and fails on divergence.
Set `WES_CALLER_REF` to pin the comparison to a caller release.
