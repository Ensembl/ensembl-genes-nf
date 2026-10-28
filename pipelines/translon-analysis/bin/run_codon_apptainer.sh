#!/usr/bin/env bash
set -euo pipefail

# Load Lmod when this script is invoked from a non-login shell.
if ! type module >/dev/null 2>&1; then
    for module_init in /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
        if [[ -r "$module_init" ]]; then
            # shellcheck disable=SC1090
            source "$module_init"
            break
        fi
    done
fi

if ! type module >/dev/null 2>&1; then
    echo "ERROR: the module command is unavailable; run this from a Codon login environment." >&2
    exit 1
fi

module load apptainer/1.4.4

for executable in apptainer nextflow java; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "ERROR: $executable is unavailable after loading Apptainer; load Nextflow and Java before running this script." >&2
        exit 1
    fi
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
pipeline_dir="$(cd -- "$script_dir/.." && pwd)"

exec nextflow run "$pipeline_dir" \
    -c "$pipeline_dir/conf/hpc_apptainer.config" \
    -profile slurm_apptainer \
    "$@"
