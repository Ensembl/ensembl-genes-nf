process DIAMOND_BLASTP {
    tag { meta.id }
    label 'process_high'
    container 'community.wave.seqera.io/library/diamond:2.1.24--61a5af76160d103f'
    publishDir "${params.outdir}/qc/diamond", mode: 'copy', overwrite: true

    input:
        tuple val(meta), path(query_protein)
        path diamond_db

    output:
        tuple val(meta), path("*_diamond.tsv"), emit: hits
        path 'versions.yml', emit: versions

    script:
        def out = "${meta.id}_diamond.tsv"
        """
        diamond blastp \
            --query ${query_protein} \
            --db ${diamond_db} \
            --threads ${task.cpus} \
            --evalue 1e-5 \
            --max-target-seqs 1 \
            --outfmt 6 qseqid sseqid pident length qstart qend qlen qcovhsp sstart send slen evalue bitscore \
            --out ${out}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: \$(diamond version | sed 's/^diamond version //')
        END_VERSIONS
        """

    stub:
        def out = "${meta.id}_diamond.tsv"
        """
        touch ${out}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: 2.1.24
        END_VERSIONS
        """
}
