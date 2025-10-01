#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/extract_geo_metadata.py \
    "Ribo-seq OR 'Ribosome Profiling' OR 'ribosome profiling' OR 'ribo-seq'" \
    -o geo_metadata.json \
    --max-results 1000 \
     \
    --max-results 1000

cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:MULTI_SOURCE_INGESTION:EXTRACT_GEO_METADATA":
    python: $(python --version | sed 's/Python //g')
    requests: $(python -c "import requests; print(requests.__version__)")
END_VERSIONS
