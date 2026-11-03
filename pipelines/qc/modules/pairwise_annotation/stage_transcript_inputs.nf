process STAGE_TRANSCRIPT_INPUTS {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_transcripts/${meta.pair_id}/inputs" }, mode: 'copy', pattern: '*.input'

    input:
    tuple val(meta), val(source_a_ref), val(source_b_ref), val(assembly_report_ref), val(evidence_a_ref), val(evidence_b_ref), val(reference_ref)
    val cache_dir

    output:
    tuple val(meta), path('source_a.input'), path('source_b.input'), path('assembly_report.input'), path('evidence_a.input'), path('evidence_b.input'), path('reference.input'), emit: staged

    script:
    """
    set -euo pipefail
    pairwise_qc.py stage --source-a '${source_a_ref}' --source-b '${source_b_ref}' --assembly-report '${assembly_report_ref}' --evidence-a '${evidence_a_ref}' --evidence-b '${evidence_b_ref}' --reference '${reference_ref}' --cache-dir '${cache_dir}' --output-source-a source_a.input --output-source-b source_b.input --output-assembly-report assembly_report.input --output-evidence-a evidence_a.input --output-evidence-b evidence_b.input --output-reference reference.input
    """

    stub:
    """
    pairwise_qc.py stage --source-a '' --source-b '' --assembly-report '' --evidence-a '' --evidence-b '' --reference '' --cache-dir '${cache_dir}' --output-source-a source_a.input --output-source-b source_b.input --output-assembly-report assembly_report.input --output-evidence-a evidence_a.input --output-evidence-b evidence_b.input --output-reference reference.input
    """
}
