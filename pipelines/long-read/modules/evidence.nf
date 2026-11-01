process COLLECT_LONG_READ_EVIDENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy'
    tag '${cohort_id}'
    input:
    path bed12s, stageAs: 'models/*'
    path normalized_manifest
    path alignment_manifests, stageAs: 'alignment_manifests/*'
    path reference_identity, stageAs: 'reference_identity.json'
    val model_records_json
    val cohort_id
    output:
    path 'evidence', emit: evidence
    path 'versions.yml', emit: versions
    script:
    """
    cat > model_records.json <<'JSON'
${model_records_json}
JSON
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json \\
        --output-dir evidence
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence
    printf '"%s":\n    evidence_packager: pipeline\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat > model_records.json <<'JSON'
${model_records_json}
JSON
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json \\
        --output-dir evidence
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence
    printf '"%s":\n    evidence_packager: stub\n' '${task.process}' > versions.yml
    """
}
