import gzip
from typing import Dict

# Module-level variables; populated by the Snakefile before any rule is evaluated.
fastq_dict: Dict = {}
probe_dict: Dict = {}
sex_dict: Dict = {}
config: Dict = {}


def get_fastq1(wildcards):
    return fastq_dict[wildcards.sample]["fq1"]


def get_fastq2(wildcards):
    return fastq_dict[wildcards.sample]["fq2"]


def get_probe_version(wildcards):
    return probe_dict[wildcards.sample]


def get_sex(wildcards):
    return sex_dict[wildcards.sample]


def parse_fastq_header(fastq_path: str, sample_name: str) -> Dict[str, str]:
    opener = (
        gzip.open(fastq_path, "rt") if fastq_path.endswith(".gz") else open(fastq_path)
    )
    with opener as f:
        header = f.readline().strip()
    try:
        parts = header.lstrip("@").split()[0].split(":")
        rgid = f"{parts[0]}:{parts[1]}"
        platform_unit = parts[2]
    except (IndexError, AttributeError):
        raise ValueError(f"Could not parse header in {fastq_path}: {header!r}")
    probe = probe_dict[sample_name]
    library_prep = config["probe_configs"][probe]["library_prep"]
    return {
        "RGID": rgid,
        "RGPU": platform_unit,
        "RGSM": sample_name,
        "RGPL": "ILLUMINA",
        "RGLB": f"{sample_name}_{library_prep}",
    }


def get_read_group_params(wildcards) -> Dict[str, str]:
    fq1_path = fastq_dict[wildcards.sample]["fq1"]
    return parse_fastq_header(fq1_path, wildcards.sample)
