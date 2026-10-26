process REPORT_UNSUPPORTED_MODEL_INPUTS {
    tag "${meta.id}:${backend}:unsupported-input"
    label 'process_light'
    publishDir path: { "${params.outdir}/per_run/${meta.id}/qc" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    input:
    tuple val(meta), val(backend)

    output:
    path("${meta.id}_${backend}_unsupported_input.tsv"), emit: report
    path 'versions.yml', emit: versions

    script:
    def technology = meta.technology ?: 'unknown'
    def representation = meta.read_representation ?: 'unknown'
    def preprocessing = meta.preprocessing_state ?: 'unknown'
    def observation_unit = meta.observation_unit ?: 'unknown'
    """
    printf 'run_accession\\tbackend\\tstatus\\treason\\ttechnology\\tread_representation\\tpreprocessing_state\\tobservation_unit\\n${meta.id}\\t${backend}\\tSKIPPED\\tunsupported_input_representation\\t${technology}\\t${representation}\\t${preprocessing}\\t${observation_unit}\\n' > ${meta.id}_${backend}_unsupported_input.tsv
    printf '"%s":\\n    report: pipeline\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'run_accession\\tbackend\\tstatus\\treason\\ttechnology\\tread_representation\\tpreprocessing_state\\tobservation_unit\\n${meta.id}\\t${backend}\\tSKIPPED\\tunsupported_input_representation\\t${meta.technology ?: 'unknown'}\\t${meta.read_representation ?: 'unknown'}\\t${meta.preprocessing_state ?: 'unknown'}\\t${meta.observation_unit ?: 'unknown'}\\n' > ${meta.id}_${backend}_unsupported_input.tsv
    printf '"%s":\\n    report: stub\\n' '${task.process}' > versions.yml
    """
}
