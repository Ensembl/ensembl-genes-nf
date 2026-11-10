process TAMA_PREPARE_SAM {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path(alignment_sam)
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path('shard.sam'), emit: prepared
    script:
    """
    python3 ${projectDir}/bin/filter_sam.py --input ${alignment_sam} --output shard.sam
    """
    stub:
    """
    touch shard.sam
    """
}

process TAMA_COLLAPSE {
    label 'process_high_memory'
    // audit container 'docker.io/sguizard/tama@sha256:66a05355155f5ecef2f585736e7fe3d9d9ad38b566f1179437828ad2d7bfd123'
    container params.use_containers ? params.tool_containers.tama : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path(shard_sam)
    path reference
    val cap_mode
    val end_mode
    val run_mode
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path(shard_sam), path('tama.bed'), emit: raw
    script:
    def args = task.ext.args ?: "-x ${cap_mode} -e ${end_mode} -rm ${run_mode}"
    """
    python3 ${projectDir}/bin/run_tama.py --sam ${shard_sam} --reference ${reference} --output tama.bed --prefix tama --args '${args}'
    """
    stub:
    """
    printf 'chr1\\t0\\t100\\ttama.1\\t0\\t+\\t0\\t100\\t0\\t1\\t100,\\t0,\\n' > tama.bed
    """
}

process TAMA_BED12_TO_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(tama_bed)
    output:
    tuple val(meta), val(shard_meta), path('native.gtf'), emit: gtf
    script:
    """
    python3 ${projectDir}/bin/bed12_to_gtf.py --bed12 ${tama_bed} --gtf native.gtf --source tama --allow-empty
    """
    stub:
    """
    printf 'chr1\\ttama\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "tama.1"; transcript_id "tama.1"; exon_number "1";\\n' > native.gtf
    """
}

process CLASSIFY_NATIVE_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}:${backend}"
    input:
    tuple val(meta), val(shard_meta), val(backend), path(native_gtf)
    val empty_message
    output:
    tuple val(meta), val(shard_meta), val(backend), path('native.gtf'), path('status.tsv'), emit: native_models
    script:
    """
    python3 ${projectDir}/bin/gtf_status.py --gtf ${native_gtf} --backend '${backend}' --status status.tsv --empty-message '${empty_message}'
    """
    stub:
    """
    printf 'chr1\\t%s\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' '${backend}' > native.gtf
    printf 'backend\\tstatus\\n%s\\tSUCCESS\\n' '${backend}' > status.tsv
    """
}

process WRITE_MODEL_INPUT_CHECKSUMS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}:${backend}"
    input:
    tuple val(meta), val(shard_meta), val(backend), path(input_files)
    output:
    tuple val(meta), val(shard_meta), val(backend), path('input_checksums.tsv'), emit: checksums
    script:
    """
    sha256sum ${input_files} > input_checksums.tsv
    """
    stub:
    """
    printf 'stub\\n' > input_checksums.tsv
    """
}

process STRINGTIE2_COLLAPSE {
    label 'process_high'
    // audit container 'quay.io/biocontainers/stringtie@sha256:84207eb013dbcff90928d8ddc75902580d3885ce8edaa36b8fab7665a3dc2e3c'
    container params.use_containers ? params.tool_containers.stringtie2 : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    path reference
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path('native.gtf'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: ''
    """
    stringtie ${shard_bam} -o native.gtf ${args}
    printf '"%s":\\n    stringtie: %s\\n' '${task.process}' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tstringtie2\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "stringtie2.1";\\n' > native.gtf
    touch versions.yml
    """
}

process STRINGTIE3_COLLAPSE {
    label 'process_high'
    // audit container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34cb224a73658053a1f31'
    container params.use_containers ? params.tool_containers.stringtie3 : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    path reference
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path('native.gtf'), emit: raw
    path 'versions.yml', emit: versions
    script:
    """
    stringtie ${shard_bam} -o native.gtf
    printf '"%s":\\n    stringtie: %s\\n' '${task.process}' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tstringtie3\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "stringtie3.1";\\n' > native.gtf
    touch versions.yml
    """
}

process SAMTOOLS_VIEW_SHARD {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}:${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    output:
    tuple val(meta), val(shard_meta), path('alignment.sam'), emit: sam
    path 'versions.yml', emit: versions
    script:
    """
    samtools view -h -F 3844 ${shard_bam} > alignment.sam
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """
    stub:
    """
    touch alignment.sam
    touch versions.yml
    """
}

process BAM_TO_ALIGNMENT_GTF {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(shard_meta), path(alignment_sam)
    output:
    tuple val(meta), val(shard_meta), path('alignment.gtf'), emit: alignment_gtf
    script:
    """
    python3 ${projectDir}/bin/bam_to_alignment_gtf.py --sam ${alignment_sam} --output alignment.gtf
    """
    stub:
    """
    printf 'chr1\\tconverter\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "read.1";\\n' > alignment.gtf
    """
}

process TMERGE_COLLAPSE {
    label 'process_high'
    // audit container 'docker.io/library/perl@sha256:587c0c65838dca2fa760e0a0ab169f7260894f208e68f38535ca644f397470c7'
    container params.use_containers ? params.tool_containers.tmerge : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(shard_meta), path(alignment_gtf)
    val tmerge_args
    output:
    tuple val(meta), val(shard_meta), path(alignment_gtf), path('native.gtf'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: tmerge_args
    """
    perl ${projectDir}/bin/tmerge ${args} ${alignment_gtf} > native.gtf
    printf '"%s":\\n    tmerge: vendored-v1.0-8b4d6e7c1c94955931946081476e326c4cece161\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'chr1\\ttmerge\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "tmerge.1";\\n' > native.gtf
    touch versions.yml
    """
}
