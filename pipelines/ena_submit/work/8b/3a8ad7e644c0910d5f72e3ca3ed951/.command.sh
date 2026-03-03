#!/bin/bash -ue
set -euo pipefail
poll_webin.py [/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/21/adb4ef7e648a67d81a1cff23fcb21c/aln2_test.queue.json, /Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/30/8f314ca4efda2dde26a94bf5771f1f/prj_GCA_000001405.28_2026_02.queue.json, /Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/3b/1452c87fb54d9df0b3c82ceb75e006/aln1_test.queue.json]       --webin-user "Webin-70684"       --webin-password "EnsemblAlignments!"       --interval 20       --max-attempts 30       --out accessions.tsv
