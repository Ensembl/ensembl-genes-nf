#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/table_formatter.py \
    processed_metadata.json \
    -o biogroup_table.csv \
    --biogroups biogroups.json \
    --format csv \
    --level biogroup \
    --include-raw \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:CROSS_MODAL_ANALYSIS:FORMAT_BIOGROUP_TABLE":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
