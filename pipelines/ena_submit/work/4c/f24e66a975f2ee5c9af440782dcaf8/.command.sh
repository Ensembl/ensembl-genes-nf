#!/bin/bash -ue
set -euo pipefail
mkdir -p xml_out

cat > row.tsv <<'EOF'
file_path	file_type	study	analysis_alias	title	description	run_accessions	run_list_path	assembly_accession	sample_accession	ref_seqs	remote_name	analysis_links	analysis_attributes	analysis_type	omit_run_refs_in_test
/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/examples/aln2.cram	cram	prj_GCA_000001405.28_2026_02	aln2_test	Example: run ERR000002 aligned to GCA_000001405.28	Per‑run manifest example row	ERR000002		GCA_000001405.28	ERS000002		aln2.cram		attr_pipeline=ensembl-genes-nf; attr_pipeline_version=dev	REFERENCE_ALIGNMENT	true
EOF

generate_analysis_xml.py       --manifest-row row.tsv       --remote-path "aln2.cram"       --md5 md5.txt       --analysis-type REFERENCE_ALIGNMENT       --omit-run-refs              --outdir xml_out

cp xml_out/analysis.xml .
cp xml_out/submission.xml .
cp xml_out/webin_submission.xml .
