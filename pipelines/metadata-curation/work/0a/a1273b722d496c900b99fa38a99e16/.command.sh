#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/quality_validator.py \
    processed_metadata.json \
    -o validation_report.json \
    --biogroups biogroups.json \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:CROSS_MODAL_ANALYSIS:VALIDATE_QUALITY":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
