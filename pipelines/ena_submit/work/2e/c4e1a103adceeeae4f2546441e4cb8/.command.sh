#!/bin/bash -ue
set -euo pipefail
curl -sS -u "Webin-70684:EnsemblAlignments!"       -H 'Content-Type: application/xml'       --data-binary @webin_submission.xml       "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit/queue" > "aln2_test.queue.json"
