process WRITE_MODEL_RECORDS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag '${cohort_id}'
    input:
    val model_records_json
    val cohort_id
    output:
    path 'model_records.json', emit: records
    script:
    """
    python3 ${projectDir}/bin/write_json.py --value '${model_records_json}' --output model_records.json
    """
    stub:
    """
    python3 ${projectDir}/bin/write_json.py --value '${model_records_json}' --output model_records.json
    """
}

process JOIN_MODEL_EVIDENCE {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${cohort_id}"
    input:
    path bed12s, stageAs: 'models/*'
    path model_records
    path direct_support_dirs, stageAs: 'direct_support/*'
    path reference
    val cohort_id
    output:
    path 'joined_support', emit: support
    script:
    """
    python3 ${projectDir}/bin/build_model_support.py --models-dir models --model-records ${model_records} --direct-support-dir direct_support --reference ${reference} --output-dir joined_support
    """
    stub:
    """
    python3 ${projectDir}/bin/build_model_support.py --models-dir models --model-records ${model_records} --direct-support-dir direct_support --reference ${reference} --output-dir joined_support
    """
}

process PACKAGE_LONG_READ_EVIDENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy'
    tag '${cohort_id}'
    input:
    path bed12s, stageAs: 'models/*'
    path normalized_manifest
    path alignment_manifests, stageAs: 'alignment_manifests/*'
    path reference_identity, stageAs: 'reference_identity.json'
    path model_records
    path direct_support_dirs, stageAs: 'direct_support/*'
    path joined_support
    val cohort_id
    output:
    path 'evidence', emit: evidence
    script:
    """
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records ${model_records} \\
        --direct-support-dir direct_support \\
        --joined-support-dir ${joined_support} \\
        --output-dir evidence
    """
    stub:
    """
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records ${model_records} \\
        --direct-support-dir direct_support \\
        --joined-support-dir ${joined_support} \\
        --output-dir evidence
    """
}

process VERIFY_LONG_READ_EVIDENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy'
    tag '${cohort_id}'
    input:
    path evidence
    val cohort_id
    output:
    path 'evidence', emit: evidence
    script:
    """
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence --require-support
    """
    stub:
    """
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence --require-support
    """
}
