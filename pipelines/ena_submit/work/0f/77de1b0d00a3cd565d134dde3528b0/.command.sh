#!/bin/bash -ue
set -euo pipefail
poll_webin.py "/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/ca/7c9a4c8aa5c805a27560cb44b08b48/prj_GCA_000001405.28_2026_02.queue.json" "/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/24/350bc7322d06ef1873293548c42e02/aln1_test.queue.json" "/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/work/17/8f4b3c6bd72e1b76bfc2732ff57438/aln2_test.queue.json"       --webin-user "Webin-70684"       --webin-password "EnsemblAlignments!"       --interval 20       --max-attempts 30       --out accessions.tsv
