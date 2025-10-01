#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/extract_sra_metadata.py \
    "test" \
    -o sra_metadata.json \
    --max-results 1 \
     \
    --max-results 1000

cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:MULTI_SOURCE_INGESTION:EXTRACT_SRA_METADATA":
    python: $(python --version | sed 's/Python //g')
    requests: $(python -c "import requests; print(requests.__version__)")
END_VERSIONS
