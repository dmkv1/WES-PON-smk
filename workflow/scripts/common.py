"""Input and param getters, keyed on the sample and its alignment units.

The caller pipeline's equivalents key on (run, sample, unit); this panel is one
flat cohort with no run dimension, so the run drops out. Nothing else differs,
and nothing here parses a FASTQ: read groups are resolved and validated once at
load time by workflow/scripts/units.py.
"""

from typing import Dict, List

# Module-level variables; populated by the Snakefile before any rule is evaluated.
units = None
unit_index: Dict = {}
units_by_sample: Dict = {}
probe_dict: Dict = {}
sex_dict: Dict = {}
samples = None
config: Dict = {}


def get_fastq1(wildcards):
    return unit_index[(wildcards.sample, wildcards.unit)]["fq1"]


def get_fastq2(wildcards):
    return unit_index[(wildcards.sample, wildcards.unit)]["fq2"]


def get_units(sample) -> List[str]:
    """Ordered unit tokens for a sample (['L001', 'L002', ...] or ['u1', ...])."""
    return units_by_sample[sample]


def get_unit_bams(wildcards) -> List[str]:
    """Per-unit query-grouped BAMs gathered by MarkDuplicates."""
    return [
        f"work/bam/units/{wildcards.sample}.{unit}.qgrp.bam"
        for unit in get_units(wildcards.sample)
    ]


def get_read_group(wildcards) -> str:
    """The unit's finished bwa '-R' argument.

    Resolved and validated once at load time by units.build_units, so nothing
    here parses a FASTQ or builds a string that could reach a shell malformed.
    """
    return unit_index[(wildcards.sample, wildcards.unit)]["rg_string"]


def get_probe_version(wildcards):
    return probe_dict[wildcards.sample]


def get_sex(wildcards):
    return sex_dict[wildcards.sample]


def get_samples(probes=None, sex=None) -> List[str]:
    """Cohort members of a probe kit, optionally of one sex.

    Every PON artifact is pooled over exactly one of these groups: the Mutect2
    panel and the CNVkit bin definitions over a kit, the CNVkit reference and
    the PureCN normal database over a kit and a sex (male and female bin sets
    differ, chrX/chrY being masked differently). Kept in one place so a rule
    cannot pool over a subtly different set than the rule feeding it.
    """
    return [
        s
        for s in probe_dict
        if (probes is None or probe_dict[s] == probes)
        and (sex is None or sex_dict[s] == sex)
    ]
