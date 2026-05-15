#!/usr/bin/env bash
set -euo pipefail

nextflow run main.nf \
    -profile slurm \
    --input inputs/lepidoptera_mir12_input.tsv \
    --outdir "${MIR12_OUTDIR:-./data_mir12}" \
    --fasta_dir "${MIR12_FASTA_DIR:-./fastas}" \
    -resume
