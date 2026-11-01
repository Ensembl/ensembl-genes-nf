process AUDIT_NATIVE_SHARD_MODELS {
    label 'process_low'
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(backend_meta), path(native_product), path(status_report)
    output:
    tuple val(meta), val(backend_meta), path("model.${meta.run_id}.${backend_meta.backend}.${backend_meta.shard_id_or_merged}.gtf"), path(status_report), path("native.${meta.run_id}.${backend_meta.backend}.${backend_meta.shard_id_or_merged}.audit.tsv"), emit: audited
    path 'versions.yml', emit: versions
    script:
    """
    if [[ '${backend_meta.status}' != 'EMPTY' ]]; then test -s ${native_product}; fi
    cp ${native_product} "model.${meta.run_id}.${backend_meta.backend}.${backend_meta.shard_id_or_merged}.gtf"
    printf 'run_id\\tbackend\\tshard_id\\tstatus\\n%s\\t%s\\t%s\\t%s\\n' '${meta.run_id}' '${backend_meta.backend}' '${backend_meta.shard_id_or_merged}' '${backend_meta.status}' > 'native.${meta.run_id}.${backend_meta.backend}.${backend_meta.shard_id_or_merged}.audit.tsv'
    printf '"%s":\\n    native_audit: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${native_product} "model.${meta.run_id}.${backend_meta.backend}.${backend_meta.shard_id_or_merged}.gtf"
    printf 'run_id\\tbackend\\tshard_id\\tstatus\\n%s\\t%s\\t%s\\t%s\\n' '${meta.run_id}' '${backend_meta.backend}' '${backend_meta.shard_id_or_merged}' '${backend_meta.status}' > 'native.${meta.run_id}.${backend_meta.backend}.${backend_meta.shard_id_or_merged}.audit.tsv'
    printf '"%s":\\n    native_audit: stub\\n' '${task.process}' > versions.yml
    """
}

process EMIT_EMPTY_NATIVE_MODEL {
    label 'process_low'
    tag "${meta.run_id}:${backend_meta.backend}"
    input:
    tuple val(meta), val(backend_meta)
    output:
    tuple val(meta), val(backend_meta), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'versions.yml', emit: versions
    script:
    """
    printf '# EMPTY: no mapped reads for run %s\\n' '${meta.run_id}' > native.gtf
    printf '%s\\tEMPTY\\tempty\\n' '${backend_meta.backend}' > status.tsv
    printf '"%s":\\n    empty_model: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '# EMPTY: no mapped reads for run %s\\n' '${meta.run_id}' > native.gtf
    printf '%s\\tEMPTY\\tempty\\n' '${backend_meta.backend}' > status.tsv
    printf '"%s":\\n    empty_model: stub\\n' '${task.process}' > versions.yml
    """
}

process AUDIT_BACKEND_SHARD_PLAN {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_run", mode: 'copy', pattern: 'backend_shard_plan.tsv', saveAs: { filename -> "${task.tag}/shards/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(shard_plan)
    val backends
    path native_audits
    output:
    path 'backend_shard_plan.tsv', emit: audit
    path 'versions.yml', emit: versions
    script:
    def backend_args = backends.join(' ')
    """
    python3 ${projectDir}/bin/audit_backend_shards.py --plan ${shard_plan} --run-id '${meta.run_id}' --backends ${backend_args} --audits ${native_audits} --report backend_shard_plan.tsv
    printf '"%s":\\n    backend_plan_audit: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${shard_plan} backend_shard_plan.tsv
    printf '"%s":\\n    backend_plan_audit: stub\\n' '${task.process}' > versions.yml
    """
}

process VALIDATE_BACKEND_SHARD_PLAN {
    label 'process_low'
    input:
    path audit
    output:
    path 'backend_shard_plan.validation.tsv', emit: validation
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/validate_backend_plan.py --audit ${audit} --output backend_shard_plan.validation.tsv
    printf '"%s":\\n    backend_plan_validation: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'status\\tvalue\\nVALID\\ttrue\\n' > backend_shard_plan.validation.tsv
    printf '"%s":\\n    backend_plan_validation: stub\\n' '${task.process}' > versions.yml
    """
}

process GTF_TO_BED12 {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename == 'versions.yml') return null
        def identity = scope == 'run' ? meta.run_id : meta.cohort_id
        def published = filename == 'models.bed12' ? "${identity}.${backend}.bed12" : filename == 'model_manifest.json' ? "${identity}.${backend}.model_manifest.json" : filename
        "${scope == 'run' ? "per_run/${meta.run_id}/models" : "per_cohort/${meta.cohort_id}/models"}/${backend}/normalized/${published}"
    }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtf)
    output:
    tuple val(meta), val(backend), val(scope), path('models.bed12'), path('model_manifest.json'), emit: bed12
    path '*.normalization_report.tsv', emit: report
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/gtf_to_bed12.py --gtf ${gtf} --bed12 models.bed12 --backend '${backend}' --scope '${scope}'
    sha256sum models.bed12 > '${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.normalization_report.tsv'
    printf '{"backend":"%s","scope":"%s","source":"%s"}\\n' '${backend}' '${scope}' '${gtf}' > model_manifest.json
    printf '"%s":\\n    gtf_to_bed12: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'chr1\\t0\\t100\\t%s.%s.model_000001\\t0\\t+\\t0\\t100\\t0\\t1\\t100,\\t0,\\n' '${backend}' '${scope}' > models.bed12
    printf '{"backend":"%s","scope":"%s","stub":true}\\n' '${backend}' '${scope}' > model_manifest.json
    printf 'stub\\n' > '${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.normalization_report.tsv'
    printf '"%s":\\n    gtf_to_bed12: stub\\n' '${task.process}' > versions.yml
    """
}

process CANONICALISE_COMBINED_MODELS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(model_manifest)
    output:
    tuple val(meta), val(backend), val(scope), path('*.canonical.bed12'), path('*.canonical.manifest.json'), emit: canonical
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/canonicalise_bed12.py --bed12 ${bed12} --manifest ${model_manifest} --output-bed12 '${backend}.${scope}.canonical.bed12' --output-manifest '${backend}.${scope}.canonical.manifest.json'
    printf '"%s":\\n    canonicaliser: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${bed12} '${backend}.${scope}.canonical.bed12'; cp ${model_manifest} '${backend}.${scope}.canonical.manifest.json'
    printf '"%s":\\n    canonicaliser: stub\\n' '${task.process}' > versions.yml
    """
}

process VALIDATE_LONG_READ_MODELS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename == 'versions.yml') return null
        def identity = scope == 'run' ? meta.run_id : meta.cohort_id
        def published = filename.endsWith('.validation.tsv') ? "${identity}.${backend}.${scope}.model_validation.tsv" : filename == "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.manifest.json" ? "${identity}.${backend}.model_manifest.json" : filename == "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.bed12" ? "${identity}.${backend}.bed12" : filename
        "${scope == 'run' ? "per_run/${meta.run_id}/models" : "per_cohort/${meta.cohort_id}/models"}/${backend}/normalized/${published}"
    }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(model_manifest)
    path reference_dictionary
    output:
    tuple val(meta), val(backend), val(scope), path("model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.bed12"), path("model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.manifest.json"), path("model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validation.tsv"), emit: validated
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/validate_bed12.py --bed12 ${bed12} --report "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validation.tsv" --reference-dictionary ${reference_dictionary}
    cp ${bed12} "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.bed12"
    cp ${model_manifest} "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.manifest.json"
    printf '"%s":\\n    bed12_validator: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${bed12} "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.bed12"; cp ${model_manifest} "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validated.manifest.json"
    printf 'metric\\tvalue\\nmodels\\t1\\nerrors\\t0\\n' > "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.validation.tsv"
    printf '"%s":\\n    bed12_validator: stub\\n' '${task.process}' > versions.yml
    """
}

process AUDIT_NATIVE_MODELS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename == 'versions.yml') return null
        def identity = scope == 'run' ? meta.run_id : meta.cohort_id
        def source_name = "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.model_summary.tsv"
        def published = filename == source_name ? "${identity}.${backend}.${scope}.model_summary.tsv" : filename
        "${scope == 'run' ? "per_run/${meta.run_id}/models" : "per_cohort/${meta.cohort_id}/models"}/${backend}/normalized/${published}"
    }
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(model_manifest), path(validation_report)
    output:
    tuple val(meta), val(backend), val(scope), path('*.model_summary.tsv'), emit: summary
    path 'versions.yml', emit: versions
    script:
    def summary_name = "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.model_summary.tsv"
    """
    python3 ${projectDir}/bin/audit_model.py --backend '${backend}' --scope '${scope}' --validation-report ${validation_report} --output '${summary_name}'
    printf '"%s":\\n    native_model_audit: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    def summary_name = "model.${meta.run_id ?: meta.cohort_id ?: meta.id}.${backend}.${scope}.model_summary.tsv"
    """
    printf 'backend\\tscope\\tmodels\\n%s\\t%s\\t1\\n' '${backend}' '${scope}' > '${summary_name}'
    printf '"%s":\\n    native_model_audit: stub\\n' '${task.process}' > versions.yml
    """
}

process COMPARE_NORMALIZED_MODELS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_cohort", mode: 'copy', saveAs: { filename -> "${cohort_id}/model_comparison/${filename}" }
    tag "${cohort_id}"
    input:
    val records_json
    val cohort_id
    output:
    path 'model_comparison.tsv', emit: comparison
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/compare_models.py --records-json '${records_json}' --output model_comparison.tsv
    printf '"%s":\\n    model_comparison: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'cohort_id\\tbackend\\tscope\\tmodel_path\\tstatus\\n' > model_comparison.tsv
    printf '"%s":\\n    model_comparison: stub\\n' '${task.process}' > versions.yml
    """
}
