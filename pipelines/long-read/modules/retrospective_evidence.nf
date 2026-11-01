process PREPARE_RETROSPECTIVE_EVIDENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/inputs", mode: 'copy', saveAs: { filename -> filename == 'reference_identity.json' ? filename : null }
    tag '${cohort_id}'
    input:
    path normalized_manifest
    path reference
    val cohort_id
    output:
    path 'models', emit: models
    path 'alignment_manifests', emit: alignment_manifests
    path 'reference_identity.json', emit: reference_identity
    path 'model_records.json', emit: model_records
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/prepare_retrospective_evidence.py \\
        --manifest ${normalized_manifest} \\
        --reference ${reference} \\
        --cohort-id ${cohort_id} \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json
    printf '"%s":\\n    retrospective_preparation: python\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p models alignment_manifests
    printf 'chr1\\t0\\t100\\tretrospective.1\\t0\\t+\\n' > models/retrospective.retrospective.run.bed12
    printf '{"run_id":"run1","retrospective":true}\\n' > alignment_manifests/run1.alignment_manifest.json
    printf '{"retrospective":true}\\n' > reference_identity.json
    printf '[]\\n' > model_records.json
    printf '"%s":\\n    retrospective_preparation: stub\\n' '${task.process}' > versions.yml
    """
}

process PACKAGE_RETROSPECTIVE_EVIDENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy'
    tag '${cohort_id}'
    input:
    path models, stageAs: 'models/*'
    path alignment_manifests, stageAs: 'alignment_manifests/*'
    path normalized_manifest
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
    printf '"%s":\\n    retrospective_packager: python\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p evidence
    printf '"%s":\\n    retrospective_packager: stub\\n' '${task.process}' > versions.yml
    """
}

process VERIFY_RETROSPECTIVE_EVIDENCE {
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
    printf '"%s":\\n    retrospective_verifier: python\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '"%s":\\n    retrospective_verifier: stub\\n' '${task.process}' > versions.yml
    """
}
