process PREPARE_BACKEND_MERGE {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtfs), val(sort_mode)
    output:
    tuple val(meta), val(backend), val(scope), path(gtfs), path('merge_input.gtf'), path('merge_gtf.status.tsv'), path('merge_status.tsv'), path('merge_skipped_inputs.tsv'), emit: prepared
    script:
    """
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${gtfs} --output merge_input.gtf --backend '${backend}' --status merge_gtf.status.tsv --skipped merge_skipped_inputs.tsv --merge-status merge_status.tsv --scope '${scope}' --sort '${sort_mode}'
    """
    stub:
    """
    cat ${gtfs} > merge_input.gtf
    printf 'backend\\tstatus\\n%s\\tCOMPLETE\\n' '${backend}' > merge_gtf.status.tsv
    printf 'backend\\tscope\\tstatus\\n%s\\t%s\\tCOMPLETE\\n' '${backend}' '${scope}' > merge_status.tsv
    printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    """
}

process MERGE_DIRECT_GTF {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(merge_input), path(status), path(merge_status), path(skipped)
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), emit: merged
    script:
    """
    cp ${merge_input} merged.gtf
    """
    stub:
    """
    cp ${merge_input} merged.gtf
    """
}

process STRINGTIE2_MERGE {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/stringtie@sha256:84207eb013dbcff90928d8ddc75902580d3885ce8edaa36b8fab7665a3dc2e3c'
    container params.use_containers ? params.tool_containers.stringtie2 : null
    tag "${meta.id}:stringtie2:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtfs), path(merge_input), path(status)
    val merge_args
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), emit: merged
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    stringtie --merge ${args} -o merged.gtf ${gtfs}
    printf '"%s":\\n    stringtie: %s\\n' '${task.process}' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    cp ${merge_input} merged.gtf
    touch versions.yml
    """
}

process STRINGTIE3_MERGE {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34cb224a73658053a1f31'
    container params.use_containers ? params.tool_containers.stringtie3 : null
    tag "${meta.id}:stringtie3:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtfs), path(merge_input), path(status)
    val merge_args
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), emit: merged
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    stringtie --merge ${args} -o merged.gtf ${gtfs}
    printf '"%s":\\n    stringtie: %s\\n' '${task.process}' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    cp ${merge_input} merged.gtf
    touch versions.yml
    """
}

process TMERGE_MERGE {
    label 'process_medium'
    // audit container 'docker.io/library/perl@sha256:587c0c65838dca2fa760e0a0ab169f7260894f208e68f38535ca644f397470c7'
    container params.use_containers ? params.tool_containers.tmerge : null
    tag "${meta.id}:tmerge:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(merge_input), path(status)
    val merge_args
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), emit: merged
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    perl ${projectDir}/bin/tmerge ${args} ${merge_input} > merged.gtf
    printf '"%s":\\n    tmerge: vendored-v1.0-8b4d6e7c1c94955931946081476e326c4cece161\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${merge_input} merged.gtf
    touch versions.yml
    """
}

process EMIT_EMPTY_MERGED_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(merge_input), path(status), path(merge_status), path(skipped)
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), emit: merged
    script:
    """
    cp ${merge_input} merged.gtf
    """
    stub:
    """
    cp ${merge_input} merged.gtf
    """
}

process WRITE_MERGE_FILELIST {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtfs), path(status)
    output:
    tuple val(meta), val(backend), val(scope), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: filelist
    path 'merge_manifest.json', emit: manifest
    script:
    """
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope '${scope}' --status-file ${status}
    """
    stub:
    """
    printf 'stub\\n' > merge_filelist.txt
    printf 'stub\\n' > merge_filelist.sha256
    printf '{"backend":"%s","scope":"%s","stub":true}\\n' '${backend}' '${scope}' > merge_manifest.json
    """
}
