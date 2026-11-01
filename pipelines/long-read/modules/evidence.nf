process WRITE_MODEL_RECORDS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag '${cohort_id}'
    input:
    val model_records_json
    val cohort_id
    output:
    path 'model_records.json', emit: records
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/write_json.py --value '${model_records_json}' --output model_records.json
    printf '"%s":\\n    model_records: python\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    python3 ${projectDir}/bin/write_json.py --value '${model_records_json}' --output model_records.json
    printf '"%s":\\n    model_records: stub\\n' '${task.process}' > versions.yml
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
    val cohort_id
    output:
    path 'evidence', emit: evidence
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records ${model_records} \\
        --output-dir evidence
    printf '"%s":\\n    evidence_packager: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records ${model_records} \\
        --output-dir evidence
    printf '"%s":\\n    evidence_packager: stub\\n' '${task.process}' > versions.yml
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
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence
    printf '"%s":\\n    evidence_verifier: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '"%s":\\n    evidence_verifier: stub\\n' '${task.process}' > versions.yml
    """
}
