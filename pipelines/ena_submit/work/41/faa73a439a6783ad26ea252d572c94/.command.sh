#!/bin/bash -ue
set -euo pipefail
mkdir -p xml_out

cat > row.tsv <<'EOF'
file_path	file_type	study	analysis_alias	title	description	run_accessions	run_list_path	assembly_accession	sample_accession	ref_seqs	remote_name	analysis_links	analysis_attributes	analysis_type	omit_run_refs_in_test
/Users/jackt/projects/ensembl-genes-nf/pipelines/ena_submit/examples/aln1.bam	bam	prj_GCA_000001405.28_2026_02	aln1_test	Example: run SRR000001 aligned to GCA_000001405.28	Per‑run manifest example row	SRR000001		GCA_000001405.28	ERS000001		aln1.bam		attr_pipeline=ensembl-genes-nf; attr_pipeline_version=dev	REFERENCE_ALIGNMENT	true
EOF

generate_analysis_xml.py       --manifest-row row.tsv       --remote-path "aln1.bam"       --md5 md5.txt       --analysis-type REFERENCE_ALIGNMENT       --omit-run-refs              --outdir xml_out

cp xml_out/analysis.xml .
cp xml_out/submission.xml .
cp xml_out/webin_submission.xml .
