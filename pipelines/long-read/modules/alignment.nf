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
    path 'reference.dict.tsv', emit: dictionary
    path 'reference_identity.json', emit: identity
    path 'versions.yml', emit: versions
    script:
    """
    minimap2 -d reference.mmi ${reference}
    awk 'BEGIN { name=""; seq_length=0 } /^>/ { if (name != "") print name "\t" seq_length; name=substr(\$1,2); seq_length=0; next } { seq_length += length(\$0) } END { if (name != "") print name "\t" seq_length }' ${reference} > reference.dict.tsv
    reference_sha=\$(sha256sum ${reference} | awk '{print \$1}')
    dictionary_sha=\$(sha256sum reference.dict.tsv | awk '{print \$1}')
    printf '{"path":"%s","sha256":"%s","dictionary_sha256":"%s"}\n' '${reference}' "\${reference_sha}" "\${dictionary_sha}" > reference_identity.json
    cat > versions.yml <<VERSIONS
"${task.process}":
    minimap2: \$(minimap2 --version)
VERSIONS
    """
    stub:
    """
    printf 'stub index\\n' > reference.mmi
    printf 'chr1\\t1000\\nchr2\\t500\\n' > reference.dict.tsv
    printf '{"path":"%s","stub":true}\\n' '${reference}' > reference_identity.json
    printf '"%s":\\n    minimap2: stub\\n' '${task.process}' > versions.yml
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
    """
    minimap2 -ax ${preset} ${secondary} ${minimap_index} ${reads} > aligned.sam
    test -s aligned.sam
    cat > versions.yml <<VERSIONS
"${task.process}":
    minimap2: \$(minimap2 --version)
VERSIONS
    """
    stub:
    """
    printf '@HD\\tVN:1.6\\tSO:unsorted\\n' > aligned.sam
    printf '"%s":\\n    minimap2: stub\\n' '${task.process}' > versions.yml
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
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.bai'), emit: alignment
    path '*.alignment_manifest.json', emit: manifest
    path 'versions.yml', emit: versions
    script:
    def secondary_json = secondary_mode ? 'true' : 'false'
    """
    samtools sort -o aligned.bam ${sam}
    samtools index aligned.bam aligned.bam.bai
    samtools quickcheck -v aligned.bam
    mv aligned.bam '${meta.run_id}.sorted.bam'
    mv aligned.bam.bai '${meta.run_id}.sorted.bam.bai'
    reference_sha=\$(sha256sum ${reference} | awk '{print \$1}')
    input_sha=\$(sha256sum ${sam} | awk '{print \$1}')
    bam_sha=\$(sha256sum '${meta.run_id}.sorted.bam' | awk '{print \$1}')
    bai_sha=\$(sha256sum '${meta.run_id}.sorted.bam.bai' | awk '{print \$1}')
    printf '{"run_id":"%s","input":"%s","reference":"%s","reference_sha256":"%s","minimap2_preset":"%s","secondary_mode":%s,"input_sha256":"%s","bam_sha256":"%s","bai_sha256":"%s"}\n' '${meta.run_id}' '${sam}' '${reference}' "\${reference_sha}" '${preset}' '${secondary_json}' "\${input_sha}" "\${bam_sha}" "\${bai_sha}" > '${meta.run_id}.alignment_manifest.json'
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'BAM\\n' > '${meta.run_id}.sorted.bam'
    printf 'BAI\\n' > '${meta.run_id}.sorted.bam.bai'
    printf '{"run_id":"%s","stub":true}\\n' '${meta.run_id}' > '${meta.run_id}.alignment_manifest.json'
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
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
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.bai'), emit: alignment
    path '*.alignment_manifest.json', emit: manifest
    path 'versions.yml', emit: versions
    script:
    """
    samtools sort -o aligned.bam ${bam}
    samtools index aligned.bam aligned.bam.bai
    samtools quickcheck -v aligned.bam
    mv aligned.bam '${meta.run_id}.sorted.bam'
    mv aligned.bam.bai '${meta.run_id}.sorted.bam.bai'
    reference_sha=\$(sha256sum ${reference} | awk '{print \$1}')
    input_sha=\$(sha256sum ${bam} | awk '{print \$1}')
    sorted_bam_sha=\$(sha256sum '${meta.run_id}.sorted.bam' | awk '{print \$1}')
    sorted_bai_sha=\$(sha256sum '${meta.run_id}.sorted.bam.bai' | awk '{print \$1}')
    printf '{"run_id":"%s","input":"%s","reference":"%s","reference_sha256":"%s","input_sha256":"%s","bam_sha256":"%s","bai_sha256":"%s"}\n' '${meta.run_id}' '${bam}' '${reference}' "\${reference_sha}" "\${input_sha}" "\${sorted_bam_sha}" "\${sorted_bai_sha}" > '${meta.run_id}.alignment_manifest.json'
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'BAM\\n' > '${meta.run_id}.sorted.bam'
    printf 'BAI\\n' > '${meta.run_id}.sorted.bam.bai'
    printf '{"run_id":"%s","stub":true}\\n' '${meta.run_id}' > '${meta.run_id}.alignment_manifest.json'
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}

process ALIGNMENT_QC {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    publishDir "${params.outdir}/per_run", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "${task.tag}/qc/alignment/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('flagstat.txt'), path('idxstats.txt'), path('quickcheck.tsv'), path('alignment_summary.tsv'), emit: qc
    path 'versions.yml', emit: versions
    script:
    """
    samtools flagstat ${bam} > flagstat.txt
    samtools idxstats ${bam} > idxstats.txt
    samtools quickcheck -v ${bam} && printf 'run_id\\tstatus\\n%s\\tVALID\\n' '${meta.run_id}' > quickcheck.tsv
    printf 'run_id\\tbam\\n%s\\t%s\\n' '${meta.run_id}' '${bam}' > alignment_summary.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    printf 'stub\\n' > flagstat.txt; printf 'chr1\\t1000\\t2\\t0\\nchr2\\t500\\t0\\t0\\n' > idxstats.txt
    printf 'run_id\\tstatus\\n%s\\tVALID\\n' '${meta.run_id}' > quickcheck.tsv
    printf 'run_id\\tbam\\n%s\\taligned.bam\\n' '${meta.run_id}' > alignment_summary.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
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
    tuple val(meta), path('alignment.sam'), emit: sam
    path 'versions.yml', emit: versions
    script:
    """
    samtools view -F 3844 ${bam} > alignment.sam
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    touch alignment.sam
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}

process EXTRACT_ALIGNMENT_FEATURES {
    label 'process_medium'
    // audit container 'docker.io/library/python@sha256:17b9724493d916f56565c131f389faff2597f46dec6ce09516f95125bc9f519f'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/per_run", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "${task.tag}/transcript_features/${filename}" }
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(sam)
    output:
    tuple val(meta), path('read_features.tsv.gz'), path('five_prime_ends.bed'), path('three_prime_ends.bed'), path('junctions.bed'), path('feature_summary.tsv'), emit: features
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/bam_to_features.py --sam ${sam} --read-features read_features.tsv.gz --five-prime five_prime_ends.bed --three-prime three_prime_ends.bed --junctions junctions.bed --summary feature_summary.tsv --run-id '${meta.run_id}'
    cat > versions.yml <<VERSIONS
"${task.process}":
    feature_extractor: pipeline
VERSIONS
    """
    stub:
    """
    printf 'read_id\\tcontig\\n' | gzip -c > read_features.tsv.gz
    touch five_prime_ends.bed three_prime_ends.bed junctions.bed
    printf 'run_id\\tfive_prime_semantics\\n%s\\talignment-oriented genomic endpoint\\n' '${meta.run_id}' > feature_summary.tsv
    printf '"%s":\\n    feature_extractor: stub\\n' '${task.process}' > versions.yml
    """
}

process VALIDATE_ALIGNMENT_REFERENCE {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), path(bam), path(bai)
    path reference_dictionary
    output:
    tuple val(meta), path(bam), path(bai), emit: alignment
    path 'reference.validation.tsv', emit: validation
    path 'versions.yml', emit: versions
    script:
    """
    samtools idxstats ${bam} > bam.idxstats.tsv
    awk -F '\\t' '\$1 != "*" { print \$1 "\\t" \$2 }' '${reference_dictionary}' > reference.contigs.tsv
    awk -F '\\t' '\$1 != "*" { print \$1 "\\t" \$2 }' bam.idxstats.tsv > bam.contigs.tsv
    if diff -u reference.contigs.tsv bam.contigs.tsv > reference.dictionary.diff; then errors=0; else errors=1; fi
    printf 'metric\\tvalue\\ncontigs\\t%s\\nerrors\\t%s\\n' "\$(wc -l < bam.contigs.tsv)" "\${errors}" > reference.validation.tsv
    if [[ "\${errors}" -ne 0 ]]; then cat reference.dictionary.diff >&2; exit 1; fi
    printf '"%s":\\n    samtools: %s\\n    reference_validation: pipeline\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'metric\\tvalue\\ncontigs\\t2\\nerrors\\t0\\n' > reference.validation.tsv
    printf '"%s":\\n    reference_validation: stub\\n' '${task.process}' > versions.yml
    """
}
