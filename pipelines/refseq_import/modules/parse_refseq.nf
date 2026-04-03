// PARSE_REFSEQ
// Parse NCBI RefSeq GFF3 → Ensembl-style GFF3.
// Optionally applies a seq-region synonym mapping (RefSeq accession → chr name).

process PARSE_REFSEQ {
    label 'process_medium'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/refseq_import", mode: 'copy', pattern: "*.refseq.gff3"

    input:
    path gff_gz        // NCBI RefSeq GFF3 (.gz)
    path synonyms_tsv  // optional: refseq_accession→chr_name TSV; pass [] to skip

    output:
    path "*.refseq.gff3", emit: gff3
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def out_name     = gff_gz.name.replaceAll(/\.gff\.gz$|\.gff3\.gz$/, '') + '.refseq.gff3'
    def syn_arg      = synonyms_tsv ? "--synonyms ${synonyms_tsv}" : ''
    def patches_arg  = params.keep_patches ? '--keep_patches' : ''
    """
    parse_refseq_gff3.py \\
        --gff3 ${gff_gz} \\
        --out  ${out_name} \\
        ${syn_arg} \\
        ${patches_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = gff_gz.name.replaceAll(/\.gff\.gz$|\.gff3\.gz$/, '') + '.refseq.gff3'
    """
    printf '##gff-version 3\\n' > ${out_name}
    printf 'chr1\\tRefSeq\\tgene\\t1000\\t50000\\t.\\t+\\t.\\tID=refseq_gene_gene-BRCA1;Name=BRCA1;biotype=protein_coding\\n' \\
        >> ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
