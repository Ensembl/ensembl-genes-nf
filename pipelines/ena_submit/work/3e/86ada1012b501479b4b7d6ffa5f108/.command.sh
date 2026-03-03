#!/bin/bash -ue
set -euo pipefail
command -v lftp >/dev/null 2>&1 || { echo 'lftp is required' >&2; exit 127; }
lftp -u Webin-70684,EnsemblAlignments! webin2.ebi.ac.uk         -e 'put aln1.bam -o aln1.bam; put md5.txt -o aln1.bam.md5; bye'
