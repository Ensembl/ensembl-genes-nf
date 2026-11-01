process INSPECT_BAM_WORKLOAD {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_run", mode: 'copy', pattern: 'workload.tsv', saveAs: { filename -> "${task.tag}/shards/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai), path(idxstats)
    output:
    tuple val(meta), path(bam), path(bai), path('workload.tsv'), emit: workload
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/inspect_bam_workload.py --idxstats ${idxstats} --output workload.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    workload_inspector: python
VERSIONS
    """
    stub:
    """
    printf 'contig\\treference_length\\tmapped_reads\\teligible\\tresource_class\\nchr1\\t1000\\t2\\ttrue\\tunassigned\\nchr2\\t500\\t0\\tfalse\\tunassigned\\n' > workload.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}

process PLAN_CONTIG_SHARDS {
    label 'process_low'
    publishDir "${params.outdir}/per_run", mode: 'copy', pattern: '*.shard_plan.tsv', saveAs: { filename -> "${task.tag}/shards/${filename}" }
    publishDir "${params.outdir}/per_run", mode: 'copy', pattern: 'shard_status.tsv', saveAs: { filename -> "${task.tag}/shards/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai), path(workload)
    val shard_mode
    val large_threshold
    val very_large_threshold
    output:
    tuple val(meta), path(bam), path(bai), path('*.shard_plan.tsv'), path('shard_plan_summary.tsv'), path('shard_status.tsv'), emit: plan
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/plan_contig_shards.py --workload ${workload} --mode ${shard_mode} --threshold ${large_threshold} --very-large-threshold ${very_large_threshold} --run-id '${meta.run_id}' --plan '${meta.run_id}.shard_plan.tsv' --summary shard_plan_summary.tsv --status shard_status.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    shard_planner: pipeline
VERSIONS
    """
    stub:
    """
    printf 'shard_id\\tcontig\\treference_length\\tmapped_reads\\teligible\\tresource_class\\tsource_run_id\\tshard_status\\n%s__chr1\\tchr1\\t1000\\t2\\ttrue\\tsmall\\t%s\\tPLANNED\\n%s__chr2\\tchr2\\t500\\t0\\tfalse\\tempty\\t%s\\tEMPTY\\n' '${meta.run_id}' '${meta.run_id}' '${meta.run_id}' '${meta.run_id}' > '${meta.run_id}.shard_plan.tsv'
    printf 'run_id\\teligible_shards\\n%s\\t1\\n' '${meta.run_id}' > shard_plan_summary.tsv
    printf 'shard_id\\tcontig\\tstatus\\n%s__chr1\\tchr1\\tPLANNED\\n%s__chr2\\tchr2\\tEMPTY\\n' '${meta.run_id}' '${meta.run_id}' > shard_status.tsv
    printf '"%s":\\n    shard_planner: stub\\n' '${task.process}' > versions.yml
    """
}

process SPLIT_BAM_BY_CONTIG {
    label 'process_high'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(shard_meta), path(bam), path(bai)
    output:
    tuple val(meta), val(shard_meta), path('*.shard.bam'), emit: shards
    path 'shard_status.tsv', emit: status
    path 'versions.yml', emit: versions
    script:
    """
    shard_prefix='${meta.run_id}.${shard_meta.shard_id}'
    if [[ '${shard_meta.contig}' == 'whole_bam' ]]; then
        samtools view -b ${bam} > "\${shard_prefix}.shard.bam"
    else
        samtools view -b ${bam} '${shard_meta.contig}' > "\${shard_prefix}.shard.bam"
    fi
    printf 'shard_id\\tcontig\\tstatus\\n%s\\t%s\\tEMITTED\\n' '${shard_meta.shard_id}' '${shard_meta.contig}' > shard_status.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    printf 'BAM\\n' > '${meta.run_id}.${shard_meta.shard_id}.shard.bam'
    printf 'shard_id\\tcontig\\tstatus\\n%s\\t%s\\tEMITTED\\n' '${shard_meta.shard_id}' '${shard_meta.contig}' > shard_status.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}

process INDEX_CONTIG_SHARD {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}:${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam)
    output:
    tuple val(meta), val(shard_meta), path('*.shard.bam'), path('*.shard.bam.bai'), emit: indexed
    path 'versions.yml', emit: versions
    script:
    """
    cp ${shard_bam} indexed.shard.bam
    samtools index indexed.shard.bam indexed.shard.bam.bai
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    cp ${shard_bam} indexed.shard.bam
    printf 'BAI\\n' > indexed.shard.bam.bai
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}

process COUNT_SOURCE_COVERAGE {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(source_bam), path(shard_plan)
    output:
    tuple val(meta), path(shard_plan), path('source_count.tsv'), emit: source_count
    path 'versions.yml', emit: versions
    script:
    """
    samtools view -c -F 4 ${source_bam} > source_count.tsv
    printf '"%s":\\n    shard_coverage: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '2\\n' > source_count.tsv
    printf '"%s":\\n    shard_coverage: stub\\n' '${task.process}' > versions.yml
    """
}

process COUNT_CONTIG_SHARD_COVERAGE {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}:${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam)
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path('shard_count.tsv'), emit: shard_count
    path 'versions.yml', emit: versions
    script:
    """
    samtools view -c -F 4 ${shard_bam} > shard_count.tsv
    printf '"%s":\\n    shard_coverage: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '2\\n' > shard_count.tsv
    printf '"%s":\\n    shard_coverage: stub\\n' '${task.process}' > versions.yml
    """
}

process COLLECT_SHARD_COVERAGE_COUNTS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(shard_plan), path(source_count), path(shard_counts), path(shard_bams)
    output:
    tuple val(meta), path(shard_plan), path(shard_bams), path('coverage_counts.tsv'), emit: counts
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/collect_shard_coverage_counts.py --source-count ${source_count} --shard-counts ${shard_counts} --output coverage_counts.tsv
    printf '"%s":\\n    shard_coverage_counts: python\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'source_mapped_reads\\t2\\nshard_mapped_reads\\t2\\n' > coverage_counts.tsv
    printf '"%s":\\n    shard_coverage_counts: stub\\n' '${task.process}' > versions.yml
    """
}

process VALIDATE_SHARD_COVERAGE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_run", mode: 'copy', pattern: 'shard.coverage.tsv', saveAs: { filename -> "${task.tag}/shards/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(shard_plan), path(shard_bams), path(coverage_counts)
    output:
    path 'shard.coverage.tsv', emit: coverage
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/validate_shard_coverage.py --plan ${shard_plan} --counts ${coverage_counts} --shard-files ${shard_bams} --report shard.coverage.tsv
    printf '"%s":\\n    shard_coverage: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'metric\\tvalue\\nplanned_eligible_shards\\t1\\nemitted_shards\\t1\\nsource_mapped_reads\\t2\\nplanned_mapped_reads\\t2\\nshard_mapped_reads\\t2\\nerrors\\t0\\n' > shard.coverage.tsv
    printf '"%s":\\n    shard_coverage: stub\\n' '${task.process}' > versions.yml
    """
}

process VALIDATE_CONTIG_SHARDS {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path('shard.validation.tsv'), emit: validated
    path 'versions.yml', emit: versions
    script:
    """
    samtools quickcheck -v ${shard_bam}
    printf 'shard_id\\tstatus\\n%s\\tVALID\\n' '${shard_meta.shard_id}' > shard.validation.tsv
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'shard_id\\tstatus\\n%s\\tVALID\\n' '${shard_meta.shard_id}' > shard.validation.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}
