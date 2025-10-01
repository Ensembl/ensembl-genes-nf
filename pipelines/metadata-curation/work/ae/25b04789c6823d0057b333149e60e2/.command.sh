#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/table_formatter.py \
    processed_metadata.json \
    -o sample_table.csv \
    --biogroups biogroups.json \
    --format csv \
    --level sample \
     \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:CROSS_MODAL_ANALYSIS:FORMAT_SAMPLE_TABLE":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
