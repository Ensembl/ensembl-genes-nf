process INDEX_BAM {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path('*.bai'), emit: indexed
    path "${meta.safe_id}.index.versions.yml", emit: versions

    script:
    """
    echo 'Creating missing BAM index for ${meta.gca_accession}:${meta.id}' >&2
    samtools quickcheck -v ${bam}
    samtools view -H ${bam} | grep -q 'SO:coordinate' || { echo 'BAM is not coordinate sorted' >&2; exit 1; }
    samtools index -@ ${task.cpus} ${bam} ${bam}.bai
    test -s ${bam}.bai || { echo 'samtools index produced no BAI' >&2; exit 1; }
    cat <<-END_VERSIONS > ${meta.safe_id}.index.versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    : > ${bam}.bai
    printf '"%s":\n    samtools: stub\n' '${task.process}' > ${meta.safe_id}.index.versions.yml
    """
}
