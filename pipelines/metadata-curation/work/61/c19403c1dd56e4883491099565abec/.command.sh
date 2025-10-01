#!/bin/bash -euo pipefail
python /Users/jackt/ensembl-genes-nf/pipelines/metadata-curation/bin/llm_context_extractor.py \
    mapped_metadata.json \
    -o context_metadata.json \


cat <<-END_VERSIONS > versions.yml
"METADATA_CURATION:ONTOLOGY_STANDARDIZATION:EXTRACT_STUDY_CONTEXT":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
