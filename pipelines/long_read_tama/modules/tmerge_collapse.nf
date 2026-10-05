process TMERGE_COLLAPSE {
    tag "${meta.id}:${shard}:tmerge"
    label 'process_high_memory'
    container 'community.wave.seqera.io/library/pip_pyfaidx_setuptools_six_pruned:b4fcc6bbecfa2bea'

    input:
    tuple val(meta), val(shard), path(reads_gtf)

    output:
    tuple val(meta), val(shard), path('*.gtf'), emit: gtf
    tuple val(meta), val(shard), path('*.status.tsv'), emit: status
    path 'versions.yml', emit: versions

    script:
    def prefix = "${meta.id}.${shard}.tmerge"
    """
    tmerge ${params.tmerge_args} --input ${reads_gtf} --output ${prefix}.gtf
    test -s ${prefix}.gtf || { echo 'tmerge collapse produced an empty GTF' >&2; exit 1; }
    model_count=\$(grep -o 'transcript_id "[^"]*"' ${prefix}.gtf | sort -u | wc -l | tr -d ' ')
    printf 'backend\\trun_accession\\tshard\\tstatus\\tmodels\\ntmerge\\t${meta.id}\\t${shard}\\tSUCCESS\\t%s\\n' "\${model_count}" > ${prefix}.status.tsv
    printf '"%s":\\n    tmerge: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    def prefix = "${meta.id}.${shard}.tmerge"
    """
    printf 'chrStub\\ttmerge\\texon\\t1\\t4\\t.\\t+\\t.\\ttranscript_id "${meta.id}.t1";\\n' > ${prefix}.gtf
    printf 'backend\\trun_accession\\tshard\\tstatus\\tmodels\\ntmerge\\t${meta.id}\\t${shard}\\tSUCCESS\\t1\\n' > ${prefix}.status.tsv
    printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
