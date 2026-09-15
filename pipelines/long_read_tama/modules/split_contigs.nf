process SPLIT_BAM_BY_CONTIG {
    tag "${meta.id}"
    label 'process_high'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('shards'), emit: shards
    path 'contig_manifest.tsv', emit: manifest

    script:
    """
    mkdir -p shards
    samtools idxstats ${bam} | awk '\$1 != "*" && \$2 > 0 {print \$1}' > contigs.txt
    while IFS= read -r contig; do
        safe=\$(printf '%s' "\$contig" | tr -c 'A-Za-z0-9._-' '_')
        samtools view -b -o "shards/${meta.id}.\${safe}.bam" ${bam} "\$contig"
        samtools index "shards/${meta.id}.\${safe}.bam"
        printf '%s\\t%s\\n' "\$contig" "shards/${meta.id}.\${safe}.bam" >> contig_manifest.tsv
    done < contigs.txt
    """

    stub:
    """
    mkdir -p shards
    touch shards/${meta.id}.stub.bam shards/${meta.id}.stub.bam.bai
    printf 'stub\\tshards/${meta.id}.stub.bam\\n' > contig_manifest.tsv
    """
}
