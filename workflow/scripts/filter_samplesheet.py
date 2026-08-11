"""Drop named samples from a PON samplesheet.

A samplesheet has one row per FASTQ pair, so a sample that failed QC occupies
one row or several. Excluding it means dropping every row whose `sample` value
matches, not just the first.

Exclusions are given by sample name, the same string the `sample` column and
every results path use. A name that matches nothing is an error: silently
building a panel that still contains the sample you meant to remove is the
failure this guards against.

Usage:
    python workflow/scripts/filter_samplesheet.py \\
        --in pon_samplesheet.csv --out pon_samplesheet_filtered.csv \\
        --exclude S1 S2 ...
    python workflow/scripts/filter_samplesheet.py \\
        --in pon_samplesheet.csv --out pon_samplesheet_filtered.csv \\
        --exclude-file docs/excluded_samples.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def read_exclusions(paths: list[Path]) -> list[str]:
    """One sample name per line. Blank lines and `#` comments are skipped."""
    names: list[str] = []
    for path in paths:
        for line in path.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                names.append(line)
    return names


def filter_samplesheet(sheet: pd.DataFrame, exclude: list[str]) -> pd.DataFrame:
    """Return `sheet` without the rows of any sample in `exclude`.

    Raises if a name matches no row, or if the result is empty.
    """
    if "sample" not in sheet.columns:
        raise ValueError("samplesheet has no 'sample' column")

    present = set(sheet["sample"])
    missing = [name for name in exclude if name not in present]
    if missing:
        raise ValueError(
            "not in the samplesheet: " + ", ".join(sorted(missing))
        )

    kept = sheet[~sheet["sample"].isin(exclude)]
    if kept.empty:
        raise ValueError("every row was excluded")
    return kept


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="infile", type=Path, required=True,
                        help="samplesheet to read")
    parser.add_argument("--out", dest="outfile", type=Path, required=True,
                        help="samplesheet to write")
    parser.add_argument("--exclude", nargs="*", default=[], metavar="SAMPLE",
                        help="sample names to drop")
    parser.add_argument("--exclude-file", nargs="*", default=[], type=Path,
                        metavar="FILE",
                        help="file of sample names to drop, one per line")
    args = parser.parse_args(argv)

    exclude = list(args.exclude) + read_exclusions(args.exclude_file)
    if not exclude:
        parser.error("nothing to exclude; give --exclude or --exclude-file")

    sheet = pd.read_csv(args.infile)
    kept = filter_samplesheet(sheet, exclude)

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    kept.to_csv(args.outfile, index=False)

    n_samples_in = sheet["sample"].nunique()
    n_samples_out = kept["sample"].nunique()
    print(
        f"{args.infile} -> {args.outfile}\n"
        f"  samples {n_samples_in} -> {n_samples_out} "
        f"({n_samples_in - n_samples_out} dropped)\n"
        f"  rows    {len(sheet)} -> {len(kept)}",
        file=sys.stderr,
    )
    for name in sorted(set(exclude)):
        print(f"  dropped {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
