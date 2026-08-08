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
