process STRINGTIE2_COLLAPSE {
    tag "${meta.id}:${shard}:stringtie2"
    label 'process_high_memory'
    conda 'bioconda::stringtie=2.2.3'
    // Use the Bioconda/Depot SIF directly: the previously pinned Wave image
    // resolves to linux/arm64 and cannot be converted on the amd64 HPC.
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0'

    input:
    tuple val(meta), val(shard), val(resource_class), val(mapped_reads), path(bam), path(bai)

    output:
    tuple val(meta), val(shard), path('*.gtf'), emit: gtf
    tuple val(meta), val(shard), path('*.status.tsv'), emit: status
    path 'versions.yml', emit: versions

    script:
    def prefix = "${meta.id}.${shard}.stringtie2"
    """
    stringtie -L -p ${task.cpus} -o ${prefix}.gtf ${bam}
    test -s ${prefix}.gtf || { echo 'StringTie2 produced an empty GTF' >&2; exit 1; }
    model_count=\$(grep -o 'transcript_id "[^"]*"' ${prefix}.gtf | sort -u | wc -l | tr -d ' ')
    printf 'backend\\trun_accession\\tshard\\tstatus\\tmodels\\nstringtie2\\t${meta.id}\\t${shard}\\tSUCCESS\\t%s\\n' "\${model_count}" > ${prefix}.status.tsv
    printf '"%s":\\n    stringtie: 2.2.3\\n' '${task.process}' > versions.yml
    """

    stub:
    def prefix = "${meta.id}.${shard}.stringtie2"
    """
    printf 'chrStub\\tstringtie\\texon\\t1\\t4\\t.\\t+\\t.\\tgene_id "${meta.id}.g1"; transcript_id "${meta.id}.t1";\\n' > ${prefix}.gtf
    printf 'backend\\trun_accession\\tshard\\tstatus\\tmodels\\nstringtie2\\t${meta.id}\\t${shard}\\tSUCCESS\\t1\\n' > ${prefix}.status.tsv
    printf '"%s":\\n    stringtie: 2.2.3-stub\\n' '${task.process}' > versions.yml
    """
}
