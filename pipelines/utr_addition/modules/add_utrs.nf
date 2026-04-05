// ADD_UTRS
// For each acceptor transcript (coding model) find a matching donor transcript
// whose exon structure contains the acceptor's CDS, then graft the donor's UTR
// exons onto the acceptor, respecting length and minimum-size limits.
//
// Donor files are supplied in priority order; the first file is searched first
// and its match wins over any match in a later file.

process ADD_UTRS {
    label 'process_medium'

    conda "conda-forge::python=3.11 bioconda::pysam=0.22 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py311h7b4e6c6_1' :
        'biocontainers/pysam:0.22.0--py311h7b4e6c6_1' }"

    publishDir path: "${params.outdir}/utr_addition", mode: 'copy', overwrite: true,
               pattern: '*.with_utrs.gff3'

    input:
    path  consolidated_gff3
    path  donor_gff3_files
    val   max_5prime
    val   max_3prime
    val   min_utr_exon

    output:
    path "*.with_utrs.gff3", emit: gff3
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_name = consolidated_gff3.name.replaceAll(/\.gff3$/, '') + '.with_utrs.gff3'
    """
    add_utrs.py \\
        --consolidated ${consolidated_gff3} \\
        --donors ${donor_gff3_files} \\
        --out ${out_name} \\
        --max-5prime ${max_5prime} \\
        --max-3prime ${max_3prime} \\
        --min-utr-exon ${min_utr_exon} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = consolidated_gff3.name.replaceAll(/\.gff3$/, '') + '.with_utrs.gff3'
    """
    printf '##gff-version 3\\n' > ${out_name}
    printf 'chr1\\tutr_addition\\tgene\\t900\\t5100\\t.\\t+\\t.\\tID=gene_00000001;biotype=protein_coding\\n' \\
        >> ${out_name}
    printf 'chr1\\tutr_addition\\tmRNA\\t900\\t5100\\t.\\t+\\t.\\tID=tx_00000001;Parent=gene_00000001;biotype=protein_coding;utr_source=stub_donor\\n' \\
        >> ${out_name}
    printf 'chr1\\tutr_addition\\texon\\t900\\t999\\t.\\t+\\t.\\tParent=tx_00000001\\n' \\
        >> ${out_name}
    printf 'chr1\\tutr_addition\\texon\\t1000\\t2000\\t.\\t+\\t.\\tParent=tx_00000001\\n' \\
        >> ${out_name}
    printf 'chr1\\tutr_addition\\tCDS\\t1000\\t2000\\t.\\t+\\t0\\tParent=tx_00000001\\n' \\
        >> ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
