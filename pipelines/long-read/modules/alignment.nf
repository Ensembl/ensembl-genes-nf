process BUILD_MINIMAP_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}/per_cohort/${params.cohort_id}/reference", mode: 'copy', saveAs: { filename -> filename in ['reference.mmi', 'reference_identity.json'] ? filename : null }
    publishDir "${params.outdir}/inputs", mode: 'copy', saveAs: { filename -> filename == 'reference_identity.json' ? filename : null }
    // audit container 'quay.io/biocontainers/minimap2@sha256:0c397895db3b494baa4f78de7110d516a1a57707d9c7df456634220bffd965ba'
    container params.use_containers ? params.tool_containers.minimap2 : null
    tag "${reference.baseName}"
    input:
    path reference
    output:
    path 'reference.mmi', emit: index
    path 'versions.yml', emit: versions
    script:
    """
    minimap2 -d reference.mmi ${reference}
    cat > versions.yml <<VERSIONS
"${task.process}":
    minimap2: \$(minimap2 --version)
VERSIONS
    """
    stub:
    """
    printf 'stub index\\n' > reference.mmi
    touch versions.yml
    """
}

process BUILD_REFERENCE_IDENTITY {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_cohort/${params.cohort_id}/reference", mode: 'copy', saveAs: { filename -> filename == 'reference_identity.json' ? filename : null }
    publishDir "${params.outdir}/inputs", mode: 'copy', saveAs: { filename -> filename == 'reference_identity.json' ? filename : null }
    tag "${reference.baseName}"
    input:
    path reference
    output:
    path 'reference.dict.tsv', emit: dictionary
    path 'reference_identity.json', emit: identity
    script:
    """
    python3 ${projectDir}/bin/reference_identity.py --reference ${reference} --dictionary reference.dict.tsv --identity reference_identity.json
    """
    stub:
    """
    printf 'chr1\\t1000\\nchr2\\t500\\n' > reference.dict.tsv
    printf '{"path":"%s","stub":true}\\n' '${reference}' > reference_identity.json
    """
}

process ALIGN_FASTQ {
    label 'process_high'
    // audit container 'quay.io/biocontainers/minimap2@sha256:0c397895db3b494baa4f78de7110d516a1a57707d9c7df456634220bffd965ba'
    container params.use_containers ? params.tool_containers.minimap2 : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(reads)
    path reference
    path minimap_index
    val preset
    val secondary_mode
    output:
    tuple val(meta), path('aligned.sam'), emit: sam
    path 'versions.yml', emit: versions
    script:
    def secondary = secondary_mode ? '' : '--secondary=no'
    def alignment_preset = meta.minimap2_preset ?: preset
    """
    minimap2 -ax ${alignment_preset} ${secondary} ${minimap_index} ${reads} > aligned.sam
    test -s aligned.sam
    cat > versions.yml <<VERSIONS
"${task.process}":
    minimap2: \$(minimap2 --version)
VERSIONS
    """
    stub:
    """
    printf '@HD\\tVN:1.6\\tSO:unsorted\\n' > aligned.sam
    touch versions.yml
    """
}

process SORT_FASTQ_ALIGNMENT {
    label 'process_high'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    publishDir "${params.outdir}/per_run", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "${task.tag}/alignment/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(sam)
    path reference
    val preset
    val secondary_mode
    output:
    tuple val(meta), path('*.sorted.bam'), path('alignment_input'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def secondary_json = secondary_mode ? 'true' : 'false'
    """
    samtools sort -o aligned.bam ${sam}
    mv aligned.bam '${meta.run_id}.sorted.bam'
    cp ${sam} alignment_input
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'BAM\\n' > '${meta.run_id}.sorted.bam'
    printf 'CSI\\n' > '${meta.run_id}.sorted.bam.csi'
    touch alignment_input
    touch versions.yml
    """
}

process SORT_INPUT_BAM {
    label 'process_high'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    publishDir "${params.outdir}/per_run", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "${task.tag}/alignment/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam)
    path reference
    output:
    tuple val(meta), path('*.sorted.bam'), path('alignment_input'), emit: raw
    path 'versions.yml', emit: versions
    script:
    """
    samtools sort -o aligned.bam ${bam}
    mv aligned.bam '${meta.run_id}.sorted.bam'
    cp ${bam} alignment_input
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'BAM\\n' > '${meta.run_id}.sorted.bam'
    printf 'CSI\\n' > '${meta.run_id}.sorted.bam.csi'
    touch alignment_input
    touch versions.yml
    """
}

process INDEX_SORTED_ALIGNMENT {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(sorted_bam), path(alignment_input)
    output:
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.csi'), path('alignment_input'), emit: indexed
    path 'versions.yml', emit: versions
    script:
    """
    cp ${sorted_bam} indexed.sorted.bam
    if [[ "${alignment_input}" != 'alignment_input' ]]; then cp ${alignment_input} alignment_input; fi
    samtools index -c -o indexed.sorted.bam.csi indexed.sorted.bam
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    cp ${sorted_bam} indexed.sorted.bam
    if [[ "${alignment_input}" != 'alignment_input' ]]; then cp ${alignment_input} alignment_input; fi
    printf 'CSI\\n' > indexed.sorted.bam.csi
    touch versions.yml
    """
}

process VALIDATE_SORTED_ALIGNMENT {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(sorted_bam), path(sorted_bai), path(alignment_input)
    output:
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.csi'), emit: alignment
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.csi'), path('alignment_input'), emit: manifest_input
    path 'versions.yml', emit: versions
    script:
    """
    cp ${sorted_bam} validated.sorted.bam
    cp ${sorted_bai} validated.sorted.bam.csi
    if [[ "${alignment_input}" != 'alignment_input' ]]; then cp ${alignment_input} alignment_input; fi
    samtools quickcheck -v validated.sorted.bam
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    cp ${sorted_bam} validated.sorted.bam
    cp ${sorted_bai} validated.sorted.bam.csi
    if [[ "${alignment_input}" != 'alignment_input' ]]; then cp ${alignment_input} alignment_input; fi
    touch versions.yml
    """
}

process WRITE_ALIGNMENT_MANIFEST {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    input:
    tuple val(meta), path(bam), path(bai), path(input)
    path reference
    val preset
    val secondary_mode
    output:
    path '*.alignment_manifest.json', emit: manifest
    script:
    def secondary = secondary_mode ? '--secondary' : ''
    def alignment_preset = meta.minimap2_preset ?: preset
    def preset_arg = alignment_preset ? "--preset '${alignment_preset}'" : ''
    """
    python3 ${projectDir}/bin/write_alignment_manifest.py --run-id '${meta.run_id}' --input ${input} --reference ${reference} --bam ${bam} --bai ${bai} ${preset_arg} ${secondary} --output '${meta.run_id}.alignment_manifest.json'
    """
    stub:
    """
    printf '{"run_id":"%s","stub":true}\\n' '${meta.run_id}' > '${meta.run_id}.alignment_manifest.json'
    """
}

process SAMTOOLS_FLAGSTAT {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    publishDir "${params.outdir}/per_run", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "${task.tag}/qc/alignment/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('flagstat.txt'), emit: flagstat
    path 'versions.yml', emit: versions
    script:
    """
    samtools flagstat ${bam} > flagstat.txt
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    printf 'stub\\n' > flagstat.txt
    touch versions.yml
    """
}

process SAMTOOLS_QC_SUMMARY {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai), path(idxstats)
    output:
    tuple val(meta), path('quickcheck.tsv'), path('alignment_summary.tsv'), emit: qc
    path 'versions.yml', emit: versions
    script:
    """
    samtools quickcheck -v ${bam}
    printf 'run_id\\tstatus\\n%s\\tVALID\\n' '${meta.run_id}' > quickcheck.tsv
    printf 'run_id\\tbam\\n%s\\t%s\\n' '${meta.run_id}' '${bam}' > alignment_summary.tsv
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'run_id\\tstatus\\n%s\\tVALID\\n' '${meta.run_id}' > quickcheck.tsv
    printf 'run_id\\tbam\\n%s\\taligned.bam\\n' '${meta.run_id}' > alignment_summary.tsv
    touch versions.yml
    """
}

process SAMTOOLS_VIEW_ALIGNMENT {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('alignment.sam.gz'), emit: sam
    path 'versions.yml', emit: versions
    script:
    """
    samtools view ${bam} | gzip -c > alignment.sam.gz
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    gzip -c </dev/null > alignment.sam.gz
    touch versions.yml
    """
}

process EXTRACT_ALIGNMENT_FEATURES {
    label 'process_medium'
    // audit container 'docker.io/library/python@sha256:17b9724493d916f56565c131f389faff2597f46dec6ce09516f95125bc9f519f'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_run", mode: 'copy', saveAs: { filename ->
        filename.endsWith('.direct_support') ? null : "${task.tag}/transcript_features/${filename}"
    }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(sam)
    output:
    tuple val(meta), path('read_features.tsv.gz'), path('five_prime_ends.bed'), path('three_prime_ends.bed'), path('junctions.bed'), path('feature_summary.tsv'), emit: features
    tuple val(meta), path("${meta.run_id}.direct_support"), emit: direct_support
    script:
    def support_arg = params.evidence_mode == 'enabled' ? '' : '--skip-support'
    """
    python3 ${projectDir}/bin/bam_to_features.py --sam ${sam} --read-features read_features.tsv.gz --five-prime five_prime_ends.bed --three-prime three_prime_ends.bed --junctions junctions.bed --summary feature_summary.tsv --support-dir '${meta.run_id}.direct_support' --run-id '${meta.run_id}' --assembly-id '${meta.assembly_id}' --sample-id '${meta.sample_id}' --library-id '${meta.library_id}' --read-set-id '${meta.read_set_id}' --technology '${meta.technology}' --min-mapq ${params.evidence_min_mapq} ${support_arg}
    """
    stub:
    def support_arg = params.evidence_mode == 'enabled' ? '' : '--skip-support'
    """
    python3 ${projectDir}/bin/bam_to_features.py --sam ${sam} --read-features read_features.tsv.gz --five-prime five_prime_ends.bed --three-prime three_prime_ends.bed --junctions junctions.bed --summary feature_summary.tsv --support-dir '${meta.run_id}.direct_support' --run-id '${meta.run_id}' --assembly-id '${meta.assembly_id}' --sample-id '${meta.sample_id}' --library-id '${meta.library_id}' --read-set-id '${meta.read_set_id}' --technology '${meta.technology}' --min-mapq ${params.evidence_min_mapq} ${support_arg}
    """
}

process SAMTOOLS_IDXSTATS {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path(bam), path(bai), path('bam.idxstats.tsv'), emit: idxstats
    path 'versions.yml', emit: versions
    script:
    """
    samtools idxstats ${bam} > bam.idxstats.tsv
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\t1000\\t2\\t0\\nchr2\\t500\\t0\\t0\\n*\\t0\\t0\\t0\\n' > bam.idxstats.tsv
    touch versions.yml
    """
}

process VALIDATE_ALIGNMENT_REFERENCE {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai), path(idxstats)
    path reference_dictionary
    output:
    tuple val(meta), path(bam), path(bai), emit: alignment
    path 'reference.validation.tsv', emit: validation
    script:
    """
    python3 ${projectDir}/bin/validate_bam_reference.py --idxstats ${idxstats} --reference-dictionary ${reference_dictionary} --report reference.validation.tsv
    """
    stub:
    """
    printf 'metric\\tvalue\\ncontigs\\t2\\nerrors\\t0\\n' > reference.validation.tsv
    """
}
