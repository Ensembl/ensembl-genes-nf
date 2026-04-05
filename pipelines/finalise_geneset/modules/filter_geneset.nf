// FILTER_GENESET
// Remove poor-quality transcripts:
//   - Single-exon models sitting within introns of multi-exon genes
//   - Transcripts with any intron shorter than min_intron_size (frameshifted)
//   - Transcripts whose CDS encodes fewer than min_orf_aa amino acids
//   - Transcripts with no supporting evidence
// If all transcripts of a gene are removed, the gene is also removed.

process FILTER_GENESET {
    label 'process_medium'

    conda "conda-forge::python=3.11 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.10.0--py311h7b4e6c6_1' :
        'biocontainers/pybedtools:0.10.0--py311h7b4e6c6_1' }"

    publishDir path: "${params.outdir}/finalise_geneset", mode: 'copy', overwrite: true,
               pattern: '*.filtered.gff3'

    input:
    path  input_gff3
    val   min_orf_aa
    val   min_intron_size

    output:
    path "*.filtered.gff3", emit: gff3
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.filtered.gff3'
    """
    filter_geneset.py \\
        --gff3 ${input_gff3} \\
        --out ${out_name} \\
        --min-orf-aa ${min_orf_aa} \\
        --min-intron-size ${min_intron_size} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.filtered.gff3'
    """
    printf '##gff-version 3\\n' > ${out_name}
    printf 'chr1\\tfinalise_geneset\\tgene\\t1000\\t5000\\t.\\t+\\t.\\tID=gene_00000001;biotype=protein_coding\\n' \\
        >> ${out_name}
    printf 'chr1\\tfinalise_geneset\\tmRNA\\t1000\\t5000\\t.\\t+\\t.\\tID=tx_00000001;Parent=gene_00000001;biotype=protein_coding\\n' \\
        >> ${out_name}
    printf 'chr1\\tfinalise_geneset\\texon\\t1000\\t2000\\t.\\t+\\t.\\tParent=tx_00000001\\n' \\
        >> ${out_name}
    printf 'chr1\\tfinalise_geneset\\texon\\t3000\\t5000\\t.\\t+\\t.\\tParent=tx_00000001\\n' \\
        >> ${out_name}
    printf 'chr1\\tfinalise_geneset\\tCDS\\t1000\\t2000\\t.\\t+\\t0\\tParent=tx_00000001\\n' \\
        >> ${out_name}
    printf 'chr1\\tfinalise_geneset\\tCDS\\t3000\\t5000\\t.\\t+\\t0\\tParent=tx_00000001\\n' \\
        >> ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
