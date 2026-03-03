#!/bin/bash -ue
set -euo pipefail
command -v lftp >/dev/null 2>&1 || { echo 'lftp is required' >&2; exit 127; }
lftp -u Webin-70684,EnsemblAlignments! webin2.ebi.ac.uk         -e 'put aln2.cram -o aln2.cram; put md5.txt -o aln2.cram.md5; bye'
