#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/biogroup_detector.py \
    processed_metadata.json \
    -o biogroups.json \
    --report biogroup_report.json \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:CROSS_MODAL_ANALYSIS:DETECT_BIOGROUPS":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
