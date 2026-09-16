process TAMA_COLLAPSE {
    tag "${meta.id}${shard ? ':' + shard : ''}:${resource_class}:${mapped_reads} reads"
    label 'process_high_memory'
    // Known per-shard TAMA failures are converted into an explicit status
    // record below. Missing executables and scheduler/resource signals remain
    // real failures.
    memory { resource_class == 'very_large' ? params.tama_memory_very_large : (resource_class == 'large' ? params.tama_memory_large : params.tama_memory_small) }
    conda 'bioconda::gs-tama=1.0.3'
    container "${params.tama_container ?: (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0' :
        'quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0')}"

    input:
    tuple val(meta), val(shard), val(resource_class), val(mapped_reads), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), val(shard), path('*_collapsed.bed'), emit: bed, optional: true
    tuple val(meta), val(shard), path('*_read.txt'), emit: read, optional: true
    tuple val(meta), val(shard), path('*_trans_read.bed'), emit: trans_read, optional: true
    tuple val(meta), val(shard), path('tama_status.tsv'), emit: status
    tuple val(meta), val(shard), path('tama_collapse.stderr'), emit: stderr
    path 'versions.yml', emit: versions

    script:
    def prefix = "${meta.id}${shard ? '.' + shard : ''}.tama"
    def args = task.ext.args ?: "-x ${params.tama_cap_mode} -e ${params.tama_end_mode} -rm ${params.tama_run_mode}"
    """
    set +e
    tama_collapse.py -s ${bam} -b BAM -f ${reference} -p ${prefix} ${args} > tama_collapse.stdout 2> tama_collapse.stderr
    rc=\$?
    set -e
    if [ \"\$rc\" -eq 0 ] && [ -s \"${prefix}_collapsed.bed\" ]; then
        printf 'accession\\tshard\\tstatus\\texit_status\\n${meta.id}\\t${shard}\\tSUCCESS\\t0\\n' > tama_status.tsv
    elif [ \"\$rc\" -eq 127 ] || [ \"\$rc\" -eq 137 ] || [ \"\$rc\" -eq 140 ] || [ \"\$rc\" -eq 143 ]; then
        cat tama_collapse.stderr >&2
        exit \"\$rc\"
    else
        rm -f \"${prefix}_collapsed.bed\"
        printf 'accession\\tshard\\tstatus\\texit_status\\n${meta.id}\\t${shard}\\tTAMA_FAILED\\t%s\\n' \"\$rc\" > tama_status.tsv
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tama_collapse: \$(tama_collapse.py -v 2>&1 | tail -n1 || true)
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}${shard ? '.' + shard : ''}.tama"
    """
    printf 'chrStub\\t0\\t4\\t${meta.id}.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${prefix}_collapsed.bed
    printf 'transcript\\tread\\n${meta.id}.1\\tstub\\n' > ${prefix}_read.txt
    printf 'chrStub\\t0\\t4\\t${meta.id}.1;stub\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${prefix}_trans_read.bed
    printf 'accession\\tshard\\tstatus\\texit_status\\n${meta.id}\\t${shard}\\tSUCCESS\\t0\\n' > tama_status.tsv
    : > tama_collapse.stderr
    printf '"%s":\\n    tama_collapse: stub\\n' '${task.process}' > versions.yml
    """
}
