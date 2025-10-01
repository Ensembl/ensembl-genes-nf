#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/evidence_metadata_processor.py \
    context_metadata.json \
    -o processed_metadata.json \
    --field-config /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/conf/field_registry.json \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:ONTOLOGY_STANDARDIZATION:PROCESS_SAMPLE_METADATA":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
