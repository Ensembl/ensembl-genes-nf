process VALIDATE_INPUT_MANIFEST {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/inputs", mode: 'copy'
    tag "${manifest.baseName}"
    input:
    path manifest
    val mode
    val cohort_id
    val manifest_root
    val ccs_mode
    val allow_unresolved_inputs
    val cache_dir
    output:
    tuple path('normalized_manifest.tsv'), path('normalized_manifest.sha256'), path('input_checksums.tsv'), path('preparation_status.tsv'), emit: normalized
    script:
    """
    python3 ${projectDir}/bin/normalize_manifest.py \
        --manifest ${manifest} --mode ${mode} --cohort-id ${cohort_id} --manifest-root ${manifest_root} --ccs-mode ${ccs_mode} ${allow_unresolved_inputs ? '--allow-unresolved-inputs' : ''} ${cache_dir ? "--cache-dir ${cache_dir}" : ''} \
        --normalized normalized_manifest.tsv --manifest-checksum normalized_manifest.sha256 \
        --checksums input_checksums.tsv --status preparation_status.tsv
    """
    stub:
    """
    python3 ${projectDir}/bin/normalize_manifest.py \\
        --manifest ${manifest} --mode ${mode} --cohort-id ${cohort_id} --manifest-root ${manifest_root} --ccs-mode ${ccs_mode} ${allow_unresolved_inputs ? '--allow-unresolved-inputs' : ''} ${cache_dir ? "--cache-dir ${cache_dir}" : ''} \\
        --normalized normalized_manifest.tsv --manifest-checksum normalized_manifest.sha256 \\
        --checksums input_checksums.tsv --status preparation_status.tsv
    """
}

process FINALIZE_PREPARED_MANIFEST {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/inputs", mode: 'copy'
    input:
    tuple path(normalized), path(manifest_checksum), path(input_checksums), path(status)
    path prepared_inputs
    output:
    tuple path('normalized_manifest.tsv'), path('normalized_manifest.sha256'), path('input_checksums.tsv'), path('preparation_status.tsv'), emit: normalized
    script:
    """
    python3 ${projectDir}/bin/finalize_prepared_manifest.py \
        --normalized ${normalized} --prepared ${prepared_inputs} \
        --output normalized_manifest.tsv --checksums input_checksums.tsv --status preparation_status.tsv
    """
    stub:
    """
    python3 ${projectDir}/bin/finalize_prepared_manifest.py \
        --normalized ${normalized} --prepared ${prepared_inputs} \
        --output normalized_manifest.tsv --checksums input_checksums.tsv --status preparation_status.tsv
    """
}

process VALIDATE_FASTQ {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path(reads), path('fastq.status.tsv'), emit: reads
    when:
    task.ext.when == null || task.ext.when
    script:
    """
    python3 ${projectDir}/bin/validate_fastq.py --reads ${reads} --run-id '${meta.id}' --status fastq.status.tsv
    """
    stub:
    """
    printf 'id\\tstatus\\n%s\\tVALID\\n' '${meta.id}' > fastq.status.tsv
    touch ${reads}
    """
}

process AUDIT_FASTQ {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(reads), path(status)
    output:
    tuple val(meta), path(reads), emit: reads
    path 'fastq.audit.tsv', emit: audit
    script:
    """
    python3 ${projectDir}/bin/audit_fastq.py --reads ${reads} --run-id '${meta.id}' --audit fastq.audit.tsv
    """
    stub:
    """
    printf 'id\\tstatus\\tbytes\\n%s\\tVALID\\t0\\n' '${meta.id}' > fastq.audit.tsv
    """
}

process VALIDATE_BAM_INPUT {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path(bam), path(bai), path('bam.status.tsv'), emit: bam
    path 'prepared_input.tsv', emit: prepared
    path 'versions.yml', emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    """
    samtools quickcheck -v ${bam}
    test -s ${bai}
    printf 'id\\tstatus\\n%s\\tVALID\\n' '${meta.id}' > bam.status.tsv
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\t%s\\t%s\\t%s\\t%s\\n' '${meta.id}' '${bam.getName()}' "\$(sha256sum ${bam} | cut -d' ' -f1)" '${bai.getName()}' "\$(sha256sum ${bai} | cut -d' ' -f1)" > prepared_input.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    printf 'id\\tstatus\\n%s\\tVALID\\n' '${meta.id}' > bam.status.tsv
    touch ${bam}
    touch ${bai}
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\t%s\\tSTUB\\t%s\\tSTUB\\n' '${meta.id}' '${bam.getName()}' '${bai.getName()}' > prepared_input.tsv
    touch versions.yml
    """
}

process STAGE_VALIDATED_INPUT {
    label 'process_light'
    tag "${meta.id}"
    input:
    tuple val(meta), path(input_artifact)
    output:
    tuple val(meta), path('staged.input'), emit: reads
    script:
    """
    cp ${input_artifact} staged.input
    """
    stub:
    """
    touch staged.input
    """
}

process STAGE_VALIDATED_FASTQ {
    label 'process_light'
    tag "${meta.id}"
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path('staged.fastq'), emit: reads
    path 'prepared_input.tsv', emit: prepared
    script:
    """
    cp ${reads} staged.fastq
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\tstaged.fastq\\t%s\\t\\t\\n' '${meta.id}' "\$(sha256sum staged.fastq | cut -d' ' -f1)" > prepared_input.tsv
    """
    stub:
    """
    touch staged.fastq
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\tstaged.fastq\\tSTUB\\t\\t\\n' '${meta.id}' > prepared_input.tsv
    """
}

process DISCOVER_LONG_READ_DATA {
    label 'process_light'
    input:
    val query
    output:
    path 'discovery.tsv', emit: discovery
    script:
    """
    printf 'query\\tstatus\\n%s\\tNOT_REQUESTED\\n' '${query}' > discovery.tsv
    """
}

process RESOLVE_LONG_READ_METADATA {
    label 'process_light'
    input:
    path manifest
    output:
    path 'resolved_metadata.tsv', emit: metadata
    script:
    """
    cp ${manifest} resolved_metadata.tsv
    """
}

process ACQUIRE_FASTQ {
    label 'process_medium'
    input:
    tuple val(meta), val(source_uri)
    output:
    tuple val(meta), path('acquired.fastq.gz'), emit: reads
    path 'acquired_input.tsv', emit: prepared
    script:
    """
    python3 ${projectDir}/bin/acquire_input.py --uri '${source_uri}' --output acquired.fastq.gz --run-id '${meta.id}' --metadata acquired_input.tsv
    """
    stub:
    """
    touch acquired.fastq.gz
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\tacquired.fastq.gz\\tSTUB\\t\\t\\n' '${meta.id}' > acquired_input.tsv
    """
}

process ACQUIRE_BAM {
    label 'process_medium'
    input:
    tuple val(meta), val(source_uri)
    output:
    tuple val(meta), path('acquired.bam'), path('acquired.bam.index'), emit: bam
    path 'acquired_input.tsv', emit: prepared
    script:
    """
    python3 ${projectDir}/bin/acquire_input.py --uri '${source_uri}' --output acquired.bam --index-uri '${meta.source_index_uri}' --index-output acquired.bam.index --run-id '${meta.id}' --metadata acquired_input.tsv
    """
    stub:
    """
    touch acquired.bam acquired.bam.index
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\tacquired.bam\\tSTUB\\tacquired.bam.index\\tSTUB\\n' '${meta.id}' > acquired_input.tsv
    """
}

process RUN_PBCCS {
    label 'process_high_memory'
    // audit container 'quay.io/biocontainers/pbccs@sha256:6876708af8f1bd3af45740d86fa2f0479bea48e72bc97118240cb9f4fbd5e13b'
    container params.use_containers ? params.tool_containers.pbccs : null
    input:
    tuple val(meta), path(subreads), path(subreads_index)
    output:
    tuple val(meta), path('ccs.bam'), emit: ccs
    path 'versions.yml', emit: versions
    script:
    """
    ccs ${subreads} ccs.bam
    cat > versions.yml <<VERSIONS
"${task.process}":
    ccs: \$(ccs --version 2>&1 | head -n1)
VERSIONS
    """
}

process INDEX_PBCCS {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/pbccs@sha256:6876708af8f1bd3af45740d86fa2f0479bea48e72bc97118240cb9f4fbd5e13b'
    container params.use_containers ? params.tool_containers.pbccs : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(ccs_bam)
    output:
    tuple val(meta), path('ccs.bam'), path('ccs.bam.pbi'), emit: ccs
    script:
    """
    pbindex ${ccs_bam}
    """
    stub:
    """
    cp ${ccs_bam} ccs.bam
    touch ccs.bam.pbi
    """
}

process BAM_TO_FASTQ {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path('converted.fastq.gz'), emit: fastq
    path 'prepared_input.tsv', emit: prepared
    path 'versions.yml', emit: versions
    script:
    """
    samtools fastq ${bam} | gzip -c > converted.fastq.gz
    test -s converted.fastq.gz
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\tconverted.fastq.gz\\t%s\\t\\t\\n' '${meta.id}' "\$(sha256sum converted.fastq.gz | cut -d' ' -f1)" > prepared_input.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    touch converted.fastq.gz
    printf 'id\\tinput\\tinput_sha256\\tinput_index\\tinput_index_sha256\\n%s\\tconverted.fastq.gz\\tSTUB\\t\\t\\n' '${meta.id}' > prepared_input.tsv
    touch versions.yml
    """
}
