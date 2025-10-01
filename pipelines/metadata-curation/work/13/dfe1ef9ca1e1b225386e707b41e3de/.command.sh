#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/ontology_mapper.py \
    geo_metadata.json \
    -o mapped_metadata.json \
    --cache-dir ontology_cache \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:ONTOLOGY_STANDARDIZATION:MAP_ONTOLOGY_TERMS":
    python: $(python --version | sed 's/Python //g')
    requests: $(python -c "import requests; print(requests.__version__)")
END_VERSIONS
