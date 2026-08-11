import pandas as pd
import os
import sys


configfile: "config.yaml"

os.environ["APPTAINER_BIND"] = config["refs"]["path"]
os.environ["SINGULARITY_BIND"] = config["refs"]["path"]

sys.path.insert(0, os.path.dirname(workflow.snakefile))
from workflow.scripts.units import (
    build_units,
    row_types_mqc_table,
    sample_renames,
    units_mqc_table,
    unit_index as _build_unit_index,
    units_by_sample as _build_units_by_sample,
)

# Resolve the samplesheet into alignment units and a validated per-sample table.
# A row naming a file pair yields one unit; a row whose fq1/fq2 are globs yields
# one per pair. Read groups are derived from the data and overridden by the
# optional flowcell/lane/library/barcode columns.
#
# units.py is vendored verbatim from WES-snakemake and must stay that way -- the
# normals in this panel are only comparable to the tumors they normalise if both
# pipelines resolve read groups identically. tests/test_vendored_units.py holds
# that line against the caller's published main.
sheet = pd.read_csv(config["samplesheet"])
units, samples, rg_warnings = build_units(
    sheet, strict=config.get("params", {}).get("rg", {}).get("strict", True)
)
# Print warnings once, not once per job.
#
# Rules with run: spawn a worker that re-imports this file, so this block runs
# again per job. 185 jobs x 167 warnings = 32550 lines of log. Ask me how I know.
#
# Env var instead of a module-level flag: the child gets a fresh namespace but
# inherits the environment, so the parent's flag survives the spawn. Fresh
# snakemake run warns again, which is what we want.
#
# Not onstart handler at it doesn't fire on --dry-run, and dry-run is when you
# actually want to read these.
if not os.environ.get("WESPON_RG_WARNED"):
    for _warning in rg_warnings:
        print(f"WARNING [read groups] {_warning}", file=sys.stderr)
    os.environ["WESPON_RG_WARNED"] = "1"

# --- Validation -----------------------------------------------------------

# The caller keys its units on (run, sample) because it groups a patient's
# samples into runs. A panel of normals is one flat cohort, so the key collapses
# to the sample -- which holds only while sample names are unique across it.
_duplicated = samples["sample"][samples["sample"].duplicated()].unique().tolist()
if _duplicated:
    raise ValueError(
        "Sample names must be unique across the cohort; this panel has no run "
        "dimension to tell two same-named samples apart:\n"
        + "\n".join(f"  {s}" for s in sorted(_duplicated))
    )

# Every capture kit in the sheet must name a probe config. The sheet carries the
# catalogue's kit token ("V6+UTR"); probe_configs is keyed by the Agilent kit
# name, which is also what the PON output paths are named after and what the
# caller's config points at. probe_configs.capture_kit joins the two.
_kit_to_probes = {
    cfg["capture_kit"]: probes for probes, cfg in config["probe_configs"].items()
}
_unknown_kits = sorted(set(samples["capture_kit"]) - set(_kit_to_probes))
if _unknown_kits:
    raise ValueError(
        f"capture_kit {', '.join(_unknown_kits)} has no probe config. Known "
        f"kits: {', '.join(sorted(_kit_to_probes))}. Add one under "
        f"probe_configs in config.yaml, with its capture_kit token."
    )

# Sex is not recoverable from a FASTQ; it has to come from the samplesheet.
# cnvkit_reference and purecn_normaldb both split on it, so an unknown one
# cannot be carried.
_bad_sex = samples[~samples["gender"].isin(["m", "f"])]
if not _bad_sex.empty:
    raise ValueError(
        "Every normal needs a known sex: the CNVkit reference and the PureCN "
        "normal database are built per sex.\n"
        + "\n".join(
            f"  {row['sample']} (patient {row['ID']}): {row['gender']!r}"
            for _, row in _bad_sex.iterrows()
        )
        + "\n\nFill the patients in catalogue/sample_annotations.tsv and "
        "regenerate the sheet."
    )

# --- Derived globals -------------------------------------------------------

outdir = config["outdir"]

unit_index = {
    (row["sample"], row["unit"]): row for row in units.to_dict("records")
}
units_by_sample = {
    sample: tokens
    for (_run, sample), tokens in _build_units_by_sample(units).items()
}

probe_dict = {row["sample"]: _kit_to_probes[row["capture_kit"]]
              for _, row in samples.iterrows()}
sex_dict = {row["sample"]: row["gender"] for _, row in samples.iterrows()}

PROBE_TYPES = sorted(set(probe_dict.values()))

# Only (probes, sex) combinations actually present in the sample sheet.
# Using the observed set avoids requesting CNVkit references for sexes that
# have no samples for a given probe type.
PROBE_SEX_COMBOS = sorted(
    {(probe_dict[s], sex_dict[s]) for s in probe_dict}
)

SAMPLES = list(probe_dict)

import workflow.scripts.common as common

common.units = units
common.unit_index = unit_index
common.units_by_sample = units_by_sample
common.probe_dict = probe_dict
common.sex_dict = sex_dict
common.samples = samples
common.config = config

from workflow.scripts.common import *


include: "workflow/rules/metadata.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/bqsr.smk"
include: "workflow/rules/coverage.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/mutect2.smk"
include: "workflow/rules/pon.smk"
include: "workflow/rules/cnvkit.smk"
include: "workflow/rules/purecn.smk"


wildcard_constraints:
    sample="[^/.]+",
    probes="[^/._]+",
    sex="[mf]",
    # Unit token: a real lane (L001) where one is known, else positional (u1).
    # '.' separates it from the sample in {sample}.{unit} paths.
    unit="L[0-9]{3}|u[0-9]+",


def _tg_notify(msg):
    env = config.get("telegram_bot_env", "")
    if not env:
        return
    run_id = config.get("run_id", "")
    prefix = f"[{run_id}] " if run_id else ""
    shell(
        f"source {env} && "
        f'curl -s -d "chat_id=$TELEGRAM_CHAT_ID&text={prefix}{msg}" '
        f'"https://api.telegram.org/bot$TELEGRAM_BOT_TOKEN/sendMessage" > /dev/null'
    )


onstart:
    _tg_notify("WES-PON-smk started 🚀")


onsuccess:
    _tg_notify("WES-PON-smk finished successfully ✅")


onerror:
    _tg_notify("WES-PON-smk FAILED ❌ — check snakemake.log")


rule all:
    input:
        # Mutect2 PON per probe type
        expand(
            f"{outdir}/PON/mutect2/{{probes}}/pon.vcf.gz",
            probes=PROBE_TYPES,
        ),
        # CNVkit PON per probe type × sex (observed combinations only)
        [
            f"{outdir}/PON/cnvkit/{probes}/reference_{sex}.cnn"
            for probes, sex in PROBE_SEX_COMBOS
        ],
        # PureCN normal database per probe type x sex (observed combinations
        # only — see PROBE_SEX_COMBOS)
        [
            f"{outdir}/PON/purecn/{probes}/normalDB_{probes}_{sex}_hg38.rds"
            for probes, sex in PROBE_SEX_COMBOS
        ],
        [
            f"{outdir}/PON/purecn/{probes}/mapping_bias_{probes}_{sex}_hg38.rds"
            for probes, sex in PROBE_SEX_COMBOS
        ],
        # Interval-consistency gate: canonical-only + identical NormalDB bins
        [
            f"{outdir}/purecn/{probes}/interval_check_{sex}.ok"
            for probes, sex in PROBE_SEX_COMBOS
        ],
        # Aggregated QC report
        f"{outdir}/qc/multiqc_report.html",
        # How each FASTQ pair's read group was decided
        f"{outdir}/metadata/units.tsv",
        f"{outdir}/metadata/samples.tsv",
        # Recalibrated BAMs
        expand(
            f"{outdir}/bam/{{sample}}/{{sample}}.bam",
            sample=SAMPLES,
        ),
