import pandas as pd
import os


configfile: "config.yaml"

os.environ["APPTAINER_BIND"] = config["refs"]["path"]
os.environ["SINGULARITY_BIND"] = config["refs"]["path"]

samples = pd.read_csv(config["samples_csv"])
samples = samples.set_index("ID", drop=False)

# --- Validation -----------------------------------------------------------

missing_files = []
for sid, row in samples.iterrows():
    for col in ["R1", "R2"]:
        if not os.path.exists(row[col]):
            missing_files.append(f"  {sid}: {row[col]}")
if missing_files:
    raise FileNotFoundError("FASTQ files not found:\n" + "\n".join(missing_files))

invalid_sex = samples[~samples["sex"].isin(["m", "f"])]
if not invalid_sex.empty:
    raise ValueError(
        "Invalid sex values in samples.csv (must be 'm' or 'f'):\n"
        + "\n".join(f"  {sid}: {row['sex']!r}" for sid, row in invalid_sex.iterrows())
    )

# --- Derived globals -------------------------------------------------------

outdir = config["outdir"]

PROBE_TYPES = samples["probes"].unique().tolist()

# Only (probes, sex) combinations actually present in the sample sheet.
# Using the observed set avoids requesting CNVkit references for sexes that
# have no samples for a given probe type.
PROBE_SEX_COMBOS = (
    samples[["probes", "sex"]]
    .drop_duplicates()
    .apply(lambda r: (r["probes"], r["sex"]), axis=1)
    .tolist()
)

fastq_dict = {
    sid: {"fq1": row["R1"], "fq2": row["R2"]} for sid, row in samples.iterrows()
}
probe_dict = {sid: row["probes"] for sid, row in samples.iterrows()}
sex_dict = {sid: row["sex"] for sid, row in samples.iterrows()}

import workflow.scripts.common as common

common.fastq_dict = fastq_dict
common.probe_dict = probe_dict
common.sex_dict = sex_dict
common.config = config

from workflow.scripts.common import *


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
        # PureCN normal database per probe type
        expand(
            f"{outdir}/PON/purecn/{{probes}}/normalDB_{{probes}}_hg38.rds",
            probes=PROBE_TYPES,
        ),
        expand(
            f"{outdir}/PON/purecn/{{probes}}/mapping_bias_{{probes}}_hg38.rds",
            probes=PROBE_TYPES,
        ),
        # Aggregated QC report
        f"{outdir}/qc/multiqc_report.html",
        # Recalibrated BAMs
        expand(
            f"{outdir}/bam/{{sample}}/{{sample}}.bam",
            sample=samples.index,
        ),
