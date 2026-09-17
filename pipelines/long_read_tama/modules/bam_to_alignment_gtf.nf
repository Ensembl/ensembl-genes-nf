process BAM_TO_ALIGNMENT_GTF {
    tag "${meta.id}:${shard}:alignment-gtf"
    label 'process_high_memory'
    // BAM conversion needs samtools; keep it out of the pip-tmerge image.
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), val(shard), val(resource_class), val(mapped_reads), path(bam), path(bai)

    output:
    tuple val(meta), val(shard), path('*.reads.gtf'), emit: gtf
    path 'versions.yml', emit: versions

    script:
    def prefix = "${meta.id}.${shard}.tmerge"
    """
    bam_to_alignment_gtf.py '${bam}' ${prefix}.reads.gtf
    test -s ${prefix}.reads.gtf || { echo 'BAM-to-GTF conversion produced no alignments' >&2; exit 1; }
    printf '"%s":\\n    samtools: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    def prefix = "${meta.id}.${shard}.tmerge"
    """
    printf 'chrStub\\ttmerge\\texon\\t1\\t4\\t.\\t+\\t.\\ttranscript_id "${meta.id}.read.1";\\n' > ${prefix}.reads.gtf
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}
