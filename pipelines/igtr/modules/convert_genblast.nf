// CONVERT_GENBLAST
// Parse GenBlast GFF output, apply PID/coverage filters, assign ig_gene/tr_gene
// biotypes from IMGT FASTA headers, and emit Ensembl-style GFF3.

process CONVERT_GENBLAST {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    tuple val(meta), path(gff)
    path  proteins    // protein batch FASTA (IMGT format headers for biotype lookup)

    output:
    tuple val(meta), path("*.filtered.gff3"), emit: gff3
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def args     = task.ext.args   ?: ''
    def min_pid  = params.min_pid      ?: 70
    def min_cov  = params.min_coverage ?: 80
    """
    convert_genblast.py \\
        --gff          ${gff} \\
        --proteins     ${proteins} \\
        --out          ${prefix}.filtered.gff3 \\
        --min-pid      ${min_pid} \\
        --min-coverage ${min_cov} \\
        --sample-id    ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    cat > ${prefix}.filtered.gff3 << 'STUB_GFF3'
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t85.5\t+\t.\tID=${prefix}_igtr_gene_00000001;Name=STUB;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=${prefix}_igtr_transcript_00000001;Parent=${prefix}_igtr_gene_00000001;Name=STUB;PID=85.50;Coverage=98.30;biotype=ig_gene
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=${prefix}_igtr_exon_00000001;Parent=${prefix}_igtr_transcript_00000001
    STUB_GFF3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
