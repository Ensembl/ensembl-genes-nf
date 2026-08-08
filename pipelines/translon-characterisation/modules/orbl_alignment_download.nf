process ORBL_ALIGNMENT_DOWNLOAD {
    label 'process_medium'
    container params.orbl_container
    tag "${meta.id}"
    publishDir "${params.outdir}/06_constraint/alignments", mode: 'copy', pattern: 'alignments/*'
    input:
    tuple val(meta), path(instances)
    output:
    tuple val(meta), path('alignments/normalised/*.fa'), emit: alignments
    tuple val(meta), path('alignments/download_report.jsonl'), emit: report
    path 'versions.yml', emit: versions
    script:
    """
    prepare_alignment_download_input.py --instances ${instances} --output alignment_requests.tsv
    mkdir -p alignments/raw
    python /opt/ORBL_tools/DownloadLocalAlignment.py ${params.orbl_alignment_set} alignment_requests.tsv alignments/raw
    normalise_orbl_alignments.py --requests alignment_requests.tsv --indir alignments/raw --outdir alignments/normalised --report alignments/download_report.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orbl_tools: \$(orbl.py --version 2>&1 | sed 's/.*Version: //')
        alignment_set: ${params.orbl_alignment_set}
    END_VERSIONS
    """
    stub:
    """
    mkdir -p alignments/normalised
    printf '>hg38 instance_id=stub\\nATGTAA\\n' > alignments/normalised/stub.fa
    printf '{"instance_id":"stub","state":"positive"}\\n' > alignments/download_report.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orbl_tools: v1.0.0
    END_VERSIONS
    """
}
