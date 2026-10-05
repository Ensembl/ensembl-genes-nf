process COMPARE_CANDIDATE_MODELS {
    tag "${meta.id}:${stage}:candidate-model-comparison"
    label 'process_light'
    publishDir "${params.outdir}/reports/comparison", mode: 'copy'

    input:
    tuple val(meta), path(model_beds), val(stage)

    output:
    path 'candidate_model_comparison_*.tsv', emit: summary
    path 'candidate_model_comparison_*.json', emit: json
    path 'candidate_model_manifest_*.tsv', emit: manifest

    script:
    def summary = "candidate_model_comparison_${stage}.tsv"
    def json = "candidate_model_comparison_${stage}.json"
    def manifest = "candidate_model_manifest_${stage}.tsv"
    """
    compare_candidate_models.py ${summary} ${json} ${model_beds.join(' ')} \\
        --manifest ${manifest} --stage '${stage}'
    """

    stub:
    def summary = "candidate_model_comparison_${stage}.tsv"
    def json = "candidate_model_comparison_${stage}.json"
    def manifest = "candidate_model_manifest_${stage}.tsv"
    """
    printf 'backend\tmodel_count\ncomparison_stub\t0\n' > ${summary}
    printf '{\"backends\": [], \"pairwise\": []}\n' > ${json}
    printf 'backend\tstage\tcanonical_model_id\tnative_model_id\n' > ${manifest}
    """
}
