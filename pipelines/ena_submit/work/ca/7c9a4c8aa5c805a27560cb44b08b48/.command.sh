#!/bin/bash -ue
set -euo pipefail
ID="prj_GCA_000001405.28_2026_02"
curl -sS -u "Webin-70684:EnsemblAlignments!"       -H 'Content-Type: application/xml'       --data-binary @webin_project.xml       "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit/queue" > "$ID.queue.json"
