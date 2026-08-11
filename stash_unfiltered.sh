#!/usr/bin/env bash
#
# Move the first PON build's pooled outputs aside so a rerun rebuilds them from
# the filtered samplesheet.
#
# Everything is moved, not deleted, into a side directory that mirrors the
# repository layout, so the old panel stays inspectable and the move is
# reversible (see --restore).
#
# Two groups move:
#
#   1. The pooled deliverables and their GenomicsDB workspaces. Under
#      --rerun-triggers mtime these are what force the pooling rules to re-run;
#      the rules re-read the samplesheet at DAG build and pick up the shorter
#      cohort. cnvkit_fix, purecn_coverage and multiqc follow on their own.
#
#   2. The QC files of the excluded samples. The multiqc rule scans results/
#      wholesale rather than its declared input list, so anything left there
#      reappears in the new report. Verified against multiqc_sources.txt: every
#      source lives under results/, nothing is parsed out of logs/, so logs are
#      left alone.
#
# Deliberately NOT moved, because they are independent of cohort membership and
# moving them would make cnvkit_autobin re-run and drag 164 cnvkit_coverage jobs
# with it: results/PON/cnvkit/*/targets.bed, antitargets.bed,
# results/cnvkit/access.bed, work/intervals/. Per-sample BAMs, VCFs, GVCFs and
# coverage files stay too -- they are still valid for the retained samples, and
# for the excluded ones they are merely unreferenced (see --with-sample-data).

set -euo pipefail
shopt -s nullglob

REPO="/mnt/data/NGS/WES/PON/WES-PON-smk"
DEST="${REPO}/results_unfiltered"
EXCLUDE_FILE="${REPO}/docs/excluded_samples.txt"
KITS=(SureSelectV6UTR SureSelectV8UTR)
SEXES=(f m)

DRY_RUN=0
WITH_SAMPLE_DATA=0
RESTORE=0

usage() {
    cat <<'EOF'
Usage: ./stash_unfiltered.sh [options]

  -n, --dry-run             print what would move, move nothing
      --with-sample-data    also stash the excluded samples' BAM/VCF/GVCF/
                            coverage outputs (large; frees disk, not required)
      --restore             move everything in the side directory back
      --dest DIR            side directory (default: results_unfiltered)
  -h, --help                this text
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--dry-run)        DRY_RUN=1 ;;
        --with-sample-data)  WITH_SAMPLE_DATA=1 ;;
        --restore)           RESTORE=1 ;;
        --dest)              DEST="$2"; shift ;;
        -h|--help)           usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

cd "$REPO"

MANIFEST="${DEST}/MANIFEST.txt"
moved=0
skipped=0

# Move one repo-relative path into DEST, keeping its directory structure.
stash() {
    local rel="$1"
    if [[ ! -e "$rel" ]]; then
        skipped=$((skipped + 1))
        return
    fi
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "  would move  $rel"
        moved=$((moved + 1))
        return
    fi
    mkdir -p "${DEST}/$(dirname "$rel")"
    if [[ -e "${DEST}/${rel}" ]]; then
        echo "  ALREADY IN DEST, leaving in place: $rel" >&2
        skipped=$((skipped + 1))
        return
    fi
    mv "$rel" "${DEST}/${rel}"
    printf '%s\n' "$rel" >> "$MANIFEST"
    echo "  moved  $rel"
    moved=$((moved + 1))
}

# ---------------------------------------------------------------- restore ----

if [[ $RESTORE -eq 1 ]]; then
    if [[ ! -f "$MANIFEST" ]]; then
        echo "no manifest at $MANIFEST -- nothing to restore" >&2
        exit 1
    fi
    while IFS= read -r rel; do
        [[ -z "$rel" ]] && continue
        if [[ ! -e "${DEST}/${rel}" ]]; then
            echo "  missing in dest, skipping: $rel" >&2
            continue
        fi
        if [[ -e "$rel" ]]; then
            echo "  exists in repo, refusing to overwrite: $rel" >&2
            continue
        fi
        if [[ $DRY_RUN -eq 1 ]]; then
            echo "  would restore  $rel"
        else
            mkdir -p "$(dirname "$rel")"
            mv "${DEST}/${rel}" "$rel"
            echo "  restored  $rel"
        fi
    done < "$MANIFEST"
    [[ $DRY_RUN -eq 0 ]] && mv "$MANIFEST" "${MANIFEST}.restored"
    echo "restore complete"
    exit 0
fi

# ------------------------------------------------------------------ stash ----

if [[ $DRY_RUN -eq 0 ]]; then
    mkdir -p "$DEST"
    {
        echo "# Stashed from ${REPO} on $(date -Is)"
        echo "# Restore with: ./stash_unfiltered.sh --restore --dest ${DEST}"
    } >> "$MANIFEST"
fi

echo "Mutect2 somatic arm"
for kit in "${KITS[@]}"; do
    stash "work/genomicsdb/${kit}"
    stash "results/PON/mutect2/${kit}/pon.vcf.gz"
    stash "results/PON/mutect2/${kit}/pon.vcf.gz.tbi"
done

echo "germline arm"
for kit in "${KITS[@]}"; do
    stash "work/genomicsdb_germline/${kit}"
    stash "results/purecn/${kit}/normals_${kit}.joint.vcf.gz"
    stash "results/purecn/${kit}/normals_${kit}.joint.vcf.gz.tbi"
done

echo "CNVkit reference"
for kit in "${KITS[@]}"; do
    for sex in "${SEXES[@]}"; do
        stash "results/PON/cnvkit/${kit}/reference_${sex}.cnn"
    done
done

echo "PureCN"
for kit in "${KITS[@]}"; do
    for sex in "${SEXES[@]}"; do
        stash "results/PON/purecn/${kit}/normalDB_${kit}_${sex}_hg38.rds"
        stash "results/PON/purecn/${kit}/mapping_bias_${kit}_${sex}_hg38.rds"
        stash "results/purecn/${kit}/coverage_files_${sex}.list"
        stash "results/purecn/${kit}/interval_check_${sex}.ok"
        # Written by the same rule as normalDB; stashed so the old and new
        # builds' diagnostics do not sit side by side under the same name.
        stash "results/purecn/${kit}/mapping_bias_hq_sites_${kit}_${sex}_hg38.bed"
        stash "results/purecn/${kit}/low_coverage_targets_${kit}_${sex}_hg38.bed"
        stash "results/purecn/${kit}/interval_weights_${kit}_${sex}_hg38.png"
    done
done

if [[ ! -f "$EXCLUDE_FILE" ]]; then
    echo "no ${EXCLUDE_FILE}; skipping the excluded samples' QC" >&2
else
    echo "QC of excluded samples"
    while IFS= read -r line; do
        sample="${line%%#*}"
        sample="$(echo "$sample" | tr -d '[:space:]')"
        [[ -z "$sample" ]] && continue
        # The trailing dot anchors the glob: P117_CTRL_BM_WES. does not match
        # P117_CTRL_BM_RESEQ_WES.
        for f in "results/qc/metrics/${sample}."*    \
                 "results/qc/fastp/${sample}."*      \
                 "results/qc/fastqc/${sample}."*; do
            stash "$f"
        done
        if [[ $WITH_SAMPLE_DATA -eq 1 ]]; then
            stash "results/bam/${sample}"
            stash "results/vcf/${sample}"
            stash "results/gvcf/${sample}"
            stash "results/coverage/${sample}"
            stash "results/purecn/coverage/${sample}.txt"
        fi
    done < "$EXCLUDE_FILE"
fi

echo
if [[ $DRY_RUN -eq 1 ]]; then
    echo "dry run: ${moved} paths would move, ${skipped} absent"
else
    echo "${moved} paths moved to ${DEST}, ${skipped} absent or already there"
    echo "manifest: ${MANIFEST}"
fi
