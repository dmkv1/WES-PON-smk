#!/usr/bin/env bash
set -euo pipefail

# Run settings (cores, memory, conda, singularity, etc.) live in a workflow
# profile under profiles/. Without one, Snakemake uses profiles/default, the
# safe floor. Select another with --workflow-profile; a --profile would be
# overridden key by key by profiles/default.
#   ./launch.sh -n                               dry run, default profile
#   ./launch.sh --workflow-profile <name>        run with profiles/<name>/
# config.yaml is loaded via the Snakefile. All args pass through to snakemake.
# config.yaml's rerun_triggers, when set, becomes --rerun-triggers unless the
# command line gives one.
trap 'rm -f snakemake.pid' EXIT

rerun_args=()
if [[ " $* " != *" --rerun-triggers"* ]]; then
    triggers=$(sed -n 's/^rerun_triggers:[[:space:]]*//p' config.yaml 2>/dev/null \
        | sed 's/[[:space:]]*#.*//; s/["'\'']//g' || true)
    read -ra triggers <<< "$triggers"
    if (( ${#triggers[@]} )); then
        rerun_args=(--rerun-triggers "${triggers[@]}")
    fi
fi

{
    snakemake "$@" "${rerun_args[@]}" &
    echo $! > snakemake.pid
    wait $!
} 2>&1 | tee snakemake.log
