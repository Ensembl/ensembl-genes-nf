process COMPARE_CANDIDATE_MODELS {
    tag "${meta.id}:${stage}:candidate-model-comparison"
    label 'process_light'
    publishDir "${params.outdir}/reports/comparison", mode: 'copy'

    input:
    // Each backend emits the same canonical BED basename. Stage the collection
    // below numbered directories so Nextflow does not reject filename collisions.
    tuple val(meta), path(model_beds, stageAs: 'bed??/*'), val(stage)

    output:
    path "candidate_model_comparison_${stage}*.tsv", emit: summary
    path "candidate_model_comparison_${stage}*.json", emit: json
    path "candidate_model_manifest_${stage}*.tsv", emit: manifest

    script:
    // Accession comparisons run once per accession; include the id in those
    // filenames so published reports cannot overwrite one another.
    def report_suffix = stage == 'accession' ? "_${meta.id}" : ''
    def summary = "candidate_model_comparison_${stage}${report_suffix}.tsv"
    def json = "candidate_model_comparison_${stage}${report_suffix}.json"
    def manifest = "candidate_model_manifest_${stage}${report_suffix}.tsv"
    """
    compare_candidate_models.py ${summary} ${json} bed??/* \\
        --manifest ${manifest} --stage '${stage}'
    """

    stub:
    def report_suffix = stage == 'accession' ? "_${meta.id}" : ''
    def summary = "candidate_model_comparison_${stage}${report_suffix}.tsv"
    def json = "candidate_model_comparison_${stage}${report_suffix}.json"
    def manifest = "candidate_model_manifest_${stage}${report_suffix}.tsv"
    """
    printf 'backend\tmodel_count\ncomparison_stub\t0\n' > ${summary}
    printf '{\"backends\": [], \"pairwise\": []}\n' > ${json}
    printf 'backend\tstage\tcanonical_model_id\tnative_model_id\n' > ${manifest}
    """
}
