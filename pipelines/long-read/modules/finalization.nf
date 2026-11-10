process FINALIZE_RUN_MANIFEST {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
    tuple path(normalized_manifest), path(input_checksums), path(preparation_status)
    path shard_plans, stageAs: 'shard_plans/*'
    path model_reports, stageAs: 'model_reports/*'
    path reference_fasta
    path reference_identity, stageAs: 'reference_identity.json'
    val run_configuration
    val stage_status
    val evidence_status
    path model_publications, stageAs: 'model_publications.jsonl'
    path merge_provenance, stageAs: 'merge_provenance.jsonl'
    output:
    path 'run_manifest.json', emit: run_manifest
    path 'final_model_checksums.tsv', emit: checksums
    script:
    """
    python3 ${projectDir}/bin/finalize_run_manifest.py \
        --configuration '${run_configuration}' --stage-status '${stage_status}' \
        --model-publications ${model_publications} --evidence-status '${evidence_status}' \
        --merge-provenance ${merge_provenance} --normalized-manifest ${normalized_manifest} \
        --input-checksums ${input_checksums} --preparation-status ${preparation_status} \
        --reference-fasta ${reference_fasta} --reference-identity reference_identity.json \
        --model-reports model_reports --run-manifest run_manifest.json --checksums final_model_checksums.tsv
    """
    stub:
    """
    printf '{"pipeline":"long-read","configuration":${run_configuration},"stage_status":${stage_status},"stub":true}\\n' > run_manifest.json
    printf 'path\\tsource_path\\tsha256\\n' > final_model_checksums.tsv
    """
}

process VERIFY_PUBLICATION_CHECKSUMS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
    tuple path(normalized_manifest), path(input_checksums), path(preparation_status)
    path shard_plans, stageAs: 'shard_plans/*'
    path model_reports, stageAs: 'model_reports/*'
    path reference_fasta
    path reference_identity, stageAs: 'reference_identity.json'
    path final_checksums
    output:
    path 'publication_checksum_verification.tsv', emit: checksum_verification
    script:
    """
    python3 ${projectDir}/bin/verify_publication_checksums.py --root . --checksums ${final_checksums} --report publication_checksum_verification.tsv
    """
    stub:
    """
    printf 'path\\tsource_path\\texpected_sha256\\tactual_sha256\\tstatus\\n# errors\\t0\\n' > publication_checksum_verification.tsv
    """
}

process COLLECT_SOFTWARE_VERSIONS {
    label 'process_light'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
    val versions_text
    output:
    path 'software_versions.yml', emit: software_versions
    script:
    def encoded_versions = versions_text.bytes.encodeBase64().toString()
    """
    python3 ${projectDir}/bin/write_text.py --base64 '${encoded_versions}' --output software_versions.yml
    """
    stub:
    """
    touch software_versions.yml
    """
}

process EXTRACT_COMBINED_TRANSCRIPTS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(model_manifest), path(validation_report)
    output:
    tuple val(meta), val(backend), val(scope), path('transcripts.fa'), emit: transcripts
    script:
    """
    python3 ${projectDir}/bin/extract_transcripts.py --bed12 ${bed12} --output transcripts.fa
    """
    stub:
    """
    printf '>stub\\nNNNN\\n' > transcripts.fa
    """
}

process PREDICT_LONGEST_ATG_ORFS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(transcripts)
    output:
    tuple val(meta), val(backend), val(scope), path('orfs.fa'), emit: orfs
    script:
    """
    python3 ${projectDir}/bin/predict_orfs.py --transcripts ${transcripts} --output orfs.fa
    """
    stub:
    """
    cp ${transcripts} orfs.fa
    """
}

process PUBLISH_ORF_PRODUCTS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "per_cohort/${meta.cohort_id}/qc/orf/${backend}/${scope}/${filename}" }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(orfs)
    output:
    path 'orfs.fa', emit: orfs
    path 'orf_manifest.tsv', emit: manifest
    script:
    """
    if [[ '${orfs}' != 'orfs.fa' ]]; then cp ${orfs} orfs.fa; fi
    printf 'backend\\tscope\\tmodel_artifact\\tproduct\\n%s\\t%s\\t%s\\torfs.fa\\n' '${backend}' '${scope}' '${orfs}' > orf_manifest.tsv
    """
    stub:
    """
    if [[ '${orfs}' != 'orfs.fa' ]]; then cp ${orfs} orfs.fa; fi
    printf 'backend\\tscope\\tmodel_artifact\\tproduct\\n%s\\t%s\\t%s\\torfs.fa\\n' '${backend}' '${scope}' '${orfs}' > orf_manifest.tsv
    """
}

process PREP_DIAMOND_DB {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/diamond@sha256:c3d56ca63d8061357ba0375d8d636519155b5601f6e93b807a827aa988f7e1d8'
    container params.use_containers ? params.tool_containers.diamond : null
    input:
    path proteins
    output:
    path 'reference.dmnd', emit: database
    path 'versions.yml', emit: versions
    script:
    """
    diamond makedb --in ${proteins} --db reference
    printf '"%s":\\n    diamond: %s\\n' '${task.process}' "\$(diamond version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'stub\\n' > reference.dmnd
    touch versions.yml
    """
}

process WRITE_DIAMOND_DATABASE_IDENTITY {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    input:
    path input_file
    val label
    output:
    path 'database_identity.json', emit: identity
    script:
    """
    python3 ${projectDir}/bin/write_file_identity.py --input ${input_file} --label '${label}' --output database_identity.json
    """
    stub:
    """
    printf '{"label":"%s","stub":true}\\n' '${label}' > database_identity.json
    """
}

process DIAMOND_BLASTP {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/diamond@sha256:c3d56ca63d8061357ba0375d8d636519155b5601f6e93b807a827aa988f7e1d8'
    container params.use_containers ? params.tool_containers.diamond : null
    input:
    tuple val(meta), val(backend), val(scope), path(orfs)
    path database
    output:
    tuple val(meta), val(backend), val(scope), path('diamond.tsv'), emit: hits
    path 'versions.yml', emit: versions
    script:
    """
    diamond blastp --query ${orfs} --db ${database.simpleName} --out diamond.tsv --outfmt 6
    printf '"%s":\\n    diamond: %s\\n' '${task.process}' "\$(diamond version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'query\\ttarget\\tscore\\n' > diamond.tsv
    touch versions.yml
    """
}

process REPORT_COMBINED_DIAMOND_MODELS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(hits)
    output:
    tuple val(meta), val(backend), val(scope), path('diamond_model_report.tsv'), emit: report
    script:
    """
    python3 ${projectDir}/bin/count_model_records.py --hits ${hits} --hits-only --report diamond_model_report.tsv
    """
    stub:
    """
    printf 'status\\tmodels\\nCOMPLETE\\t0\\n' > diamond_model_report.tsv
    """
}

process COMPARE_CANDIDATE_MODELS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "per_cohort/${meta.cohort_id}/qc/diamond/${backend}/${scope}/${filename}" }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(diamond_report)
    output:
    path 'candidate_model_comparison.tsv', emit: comparison
    script:
    """
    python3 ${projectDir}/bin/count_model_records.py --models ${bed12} --hits ${diamond_report} --backend '${backend}' --scope '${scope}' --report candidate_model_comparison.tsv
    """
    stub:
    """
    printf 'backend\\tscope\\tmodels\\tdiamond_hits\\tclassification_policy\\n%s\\t%s\\t1\\t0\\tlongest_atg_orf\\n' '${backend}' '${scope}' > candidate_model_comparison.tsv
    """
}
