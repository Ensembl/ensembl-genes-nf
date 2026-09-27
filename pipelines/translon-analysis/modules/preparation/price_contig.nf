process MAKE_PRICE_CONTIG_MANIFEST {
    tag "${meta.id}:price-contigs"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(fai), path(bam)

    output:
    tuple val(meta), path('price_contigs.tsv'), emit: manifest

    script:
    """
    awk 'BEGIN { OFS="\\t" } { print \$1, \$2 }' ${fai} > price_contigs.tsv
    test -s price_contigs.tsv
    """

    stub:
    """
    printf 'chr20\\t64444167\\n' > price_contigs.tsv
    """
}

process PREPARE_PRICE_CONTIG {
    tag "${meta.id}:${contig}"
    label 'process_high'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam), path(bai), path(gtf), path(fasta), path(fai), val(contig)

    output:
    tuple val(meta), path('price.bam'), path('price.bam.bai'), path('price.gtf'), path('price.fa'), emit: contig

    script:
    """
    samtools view -bh ${bam} '${contig}' > price.bam
    samtools index price.bam
    samtools faidx ${fasta} '${contig}' > price.fa
    awk -v c='${contig}' '($0 ~ /^#/ || $1 == c)' ${gtf} > price.gtf
    test -s price.fa && test -s price.gtf
    """

    stub:
    """
    touch price.bam price.bam.bai
    printf '>stub\\nN\\n' > price.fa
    printf '# stub\\n' > price.gtf
    """
}
