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
    val model_publications
    val merge_provenance
    output:
    path 'run_manifest.json', emit: run_manifest
    path 'final_model_checksums.tsv', emit: checksums
    path 'publication_checksum_verification.tsv', emit: checksum_verification
    script:
    """
    python3 ${projectDir}/bin/finalize_run_manifest.py \
        --configuration '${run_configuration}' --stage-status '${stage_status}' \
        --model-publications '${model_publications}' --evidence-status '${evidence_status}' \
        --merge-provenance '${merge_provenance}' --normalized-manifest ${normalized_manifest} \
        --input-checksums ${input_checksums} --preparation-status ${preparation_status} \
        --reference-fasta ${reference_fasta} --reference-identity reference_identity.json \
        --model-reports model_reports --run-manifest run_manifest.json --checksums final_model_checksums.tsv
    python3 ${projectDir}/bin/verify_publication_checksums.py --root . --checksums final_model_checksums.tsv --report publication_checksum_verification.tsv
    """
    stub:
    """
    printf '{"pipeline":"long-read","configuration":${run_configuration},"stage_status":${stage_status},"stub":true}\\n' > run_manifest.json
    printf 'path\\tsource_path\\tsha256\\n' > final_model_checksums.tsv
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
    """
    cat > versions.input.yml <<'VERSIONS'
${versions_text}
VERSIONS
    cp versions.input.yml software_versions.yml
    """
    stub:
    """
    printf 'long_read: stub\\n' > software_versions.yml
    """
}

process EXTRACT_COMBINED_TRANSCRIPTS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(model_manifest), path(validation_report)
    output:
    tuple val(meta), val(backend), val(scope), path('transcripts.fa'), emit: transcripts
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/extract_transcripts.py --bed12 ${bed12} --output transcripts.fa
    printf '"%s":\\n    transcript_extractor: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '>stub\\nNNNN\\n' > transcripts.fa
    printf '"%s":\\n    transcript_extractor: stub\\n' '${task.process}' > versions.yml
    """
}

process PREDICT_LONGEST_ATG_ORFS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(transcripts)
    output:
    tuple val(meta), val(backend), val(scope), path('orfs.fa'), emit: orfs
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/predict_orfs.py --transcripts ${transcripts} --output orfs.fa
    printf '"%s":\\n    orf_predictor: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${transcripts} orfs.fa
    printf '"%s":\\n    orf_predictor: stub\\n' '${task.process}' > versions.yml
    """
}

process PUBLISH_ORF_PRODUCTS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "per_cohort/${meta.cohort_id}/qc/orf/${backend}/${scope}/${filename}" }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(orfs)
    output:
    path 'orfs.fa', emit: orfs
    path 'orf_manifest.tsv', emit: manifest
    path 'versions.yml', emit: versions
    script:
    """
    if [[ '${orfs}' != 'orfs.fa' ]]; then cp ${orfs} orfs.fa; fi
    printf 'backend\\tscope\\tmodel_artifact\\tproduct\\n%s\\t%s\\t%s\\torfs.fa\\n' '${backend}' '${scope}' '${orfs}' > orf_manifest.tsv
    printf '"%s":\\n    orf_product: longest_atg\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    if [[ '${orfs}' != 'orfs.fa' ]]; then cp ${orfs} orfs.fa; fi
    printf 'backend\\tscope\\tmodel_artifact\\tproduct\\n%s\\t%s\\t%s\\torfs.fa\\n' '${backend}' '${scope}' '${orfs}' > orf_manifest.tsv
    printf '"%s":\\n    orf_product: stub\\n' '${task.process}' > versions.yml
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
    path 'database_identity.json', emit: identity
    path 'versions.yml', emit: versions
    script:
    """
    diamond makedb --in ${proteins} --db reference
    sha256sum ${proteins} > database_identity.json
    printf '"%s":\\n    diamond: %s\\n' '${task.process}' "\$(diamond version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'stub\\n' > reference.dmnd; printf '{"stub":true}\\n' > database_identity.json
    printf '"%s":\\n    diamond: stub\\n' '${task.process}' > versions.yml
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
    path 'database_identity.tsv', emit: database_identity
    path 'versions.yml', emit: versions
    script:
    """
    diamond blastp --query ${orfs} --db ${database.simpleName} --out diamond.tsv --outfmt 6
    sha256sum ${database} > database_identity.tsv
    printf '"%s":\\n    diamond: %s\\n' '${task.process}' "\$(diamond version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'query\\ttarget\\tscore\\n' > diamond.tsv
    printf 'stub\\n' > database_identity.tsv
    printf '"%s":\\n    diamond: stub\\n' '${task.process}' > versions.yml
    """
}

process REPORT_COMBINED_DIAMOND_MODELS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(hits)
    output:
    tuple val(meta), val(backend), val(scope), path('diamond_model_report.tsv'), emit: report
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/count_model_records.py --hits ${hits} --hits-only --report diamond_model_report.tsv
    printf '"%s":\\n    diamond_report: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'status\\tmodels\\nCOMPLETE\\t0\\n' > diamond_model_report.tsv
    printf '"%s":\\n    diamond_report: stub\\n' '${task.process}' > versions.yml
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
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/count_model_records.py --models ${bed12} --hits ${diamond_report} --backend '${backend}' --scope '${scope}' --report candidate_model_comparison.tsv
    printf '"%s":\\n    candidate_comparison: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'backend\\tscope\\tmodels\\tdiamond_hits\\tclassification_policy\\n%s\\t%s\\t1\\t0\\tlongest_atg_orf\\n' '${backend}' '${scope}' > candidate_model_comparison.tsv
    printf '"%s":\\n    candidate_comparison: stub\\n' '${task.process}' > versions.yml
    """
}
