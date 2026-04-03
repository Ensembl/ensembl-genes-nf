// FILTER_NCRNA
// Filter cmsearch tblout (Rfam) and optional BLASTN tabular (miRBase) hits.
// Assigns ncRNA biotypes and writes GFF3 with gene/transcript/exon hierarchy.

process FILTER_NCRNA {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/short_ncrna", mode: 'copy', pattern: "*.ncrna.gff3"

    input:
    tuple val(meta), path(tblout)
    path  blast_tsv   // optional; pass [] to skip miRNA BLAST

    output:
    tuple val(meta), path("*.ncrna.gff3"), emit: gff3
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: meta.id
    def args       = task.ext.args   ?: ''
    def blast_arg  = blast_tsv ? "--blast ${blast_tsv}" : ''
    def max_eval   = params.cmsearch_max_evalue   ?: 0.01
    def min_score  = params.cmsearch_min_score     ?: 0
    def blast_pid  = params.mirna_blast_min_pid    ?: 80
    def blast_eval = params.mirna_blast_max_evalue ?: 0.01
    """
    filter_ncrna.py \\
        --tblout          ${tblout} \\
        ${blast_arg} \\
        --out             ${prefix}.ncrna.gff3 \\
        --min-score       ${min_score} \\
        --max-evalue      ${max_eval} \\
        --blast-min-pid   ${blast_pid} \\
        --blast-max-evalue ${blast_eval} \\
        --sample-id       ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\nchr1\tcmsearch\tgene\t1000\t1088\t87.3\t+\t.\tID=${prefix}_ncrna_gene_00000001;Name=mir-21;biotype=miRNA\n' \\
        > ${prefix}.ncrna.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
