process DIAMOND_BLASTP {
    tag { "${meta.id}:${meta.backend ?: 'default'}:${meta.scope ?: 'cohort'}" }
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
        // Use explicit map indexing here.  In the long-read QC path `meta` is
        // a map containing both cohort and backend; interpolating the map (or
        // relying on ambiguous property resolution) produces a filename such
        // as `[id:..., backend:...]_diamond.tsv`, which Diamond parses as
        // multiple arguments to --out.
        def sample_id = meta['id'].toString()
        def backend = meta.containsKey('backend') ? meta['backend'].toString() : 'default'
        def scope = meta.containsKey('scope') ? meta['scope'].toString() : 'cohort'
        def out = "${sample_id}_${scope}_${backend}_diamond.tsv"
        """
        if grep -q '^>' ${query_protein}; then
            diamond blastp \
                --query ${query_protein} \
                --db ${diamond_db} \
                --threads ${task.cpus} \
                --evalue 1e-5 \
                --max-target-seqs 1 \
                --outfmt 6 qseqid sseqid pident length qstart qend qlen qcovhsp sstart send slen evalue bitscore \
                --out ${out}
        else
            : > ${out}
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: \$(diamond version | sed 's/^diamond version //')
        END_VERSIONS
        """

    stub:
        def sample_id = meta['id'].toString()
        def backend = meta.containsKey('backend') ? meta['backend'].toString() : 'default'
        def scope = meta.containsKey('scope') ? meta['scope'].toString() : 'cohort'
        def out = "${sample_id}_${scope}_${backend}_diamond.tsv"
        """
        touch ${out}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: 2.1.24
        END_VERSIONS
        """
}
