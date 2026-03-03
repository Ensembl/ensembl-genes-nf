#!/bin/bash -ue
set -euo pipefail
poll_webin.py "/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/30/8f314ca4efda2dde26a94bf5771f1f/prj_GCA_000001405.28_2026_02.queue.json" "/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/ea/3cfc8653137810366c512728c97828/aln1_test.queue.json" "/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/a9/fe98a4234a542725da9cfb03a82055/aln2_test.queue.json"       --webin-user "Webin-70684"       --webin-password "EnsemblAlignments!"       --interval 20       --max-attempts 30       --out accessions.tsv
