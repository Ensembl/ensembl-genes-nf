process COLLECT_RETROSPECTIVE_EVIDENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> filename == 'evidence' || filename.startsWith('evidence/') ? filename : null }
    publishDir "${params.outdir}/inputs", mode: 'copy', saveAs: { filename -> filename == 'reference_identity.json' ? filename : null }
    tag '${cohort_id}'
    input:
    path normalized_manifest
    path reference
    val cohort_id
    output:
    path 'evidence', emit: evidence
    path 'reference_identity.json', emit: reference_identity
    path 'versions.yml', emit: versions
    script:
    """
    mkdir -p models alignment_manifests
    python3 ${projectDir}/bin/prepare_retrospective_evidence.py \\
        --manifest ${normalized_manifest} \\
        --reference ${reference} \\
        --cohort-id ${cohort_id} \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json \\
        --output-dir evidence
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence
    printf 'retrospective_evidence\tpipeline\n' > versions.yml
    """
    stub:
    """
    mkdir -p models alignment_manifests
    python3 ${projectDir}/bin/prepare_retrospective_evidence.py \\
        --manifest ${normalized_manifest} \\
        --reference ${reference} \\
        --cohort-id ${cohort_id} \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json
    python3 ${projectDir}/bin/package_evidence.py \\
        --models-dir models \\
        --alignment-manifests-dir alignment_manifests \\
        --normalized-manifest ${normalized_manifest} \\
        --reference-identity reference_identity.json \\
        --model-records model_records.json \\
        --output-dir evidence
    python3 ${projectDir}/bin/verify_evidence_publication.py --evidence-dir evidence
    printf 'retrospective_evidence\tstub\n' > versions.yml
    """
}
