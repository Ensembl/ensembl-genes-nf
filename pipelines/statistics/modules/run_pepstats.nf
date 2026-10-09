process RUN_PEPSTATS {
    label 'process_low'

    tag "${meta.dbname}:pepstats"
    label 'pepstats'

    publishDir { "${params.outdir}/${meta.gca}" }, mode: 'copy'

    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path("pepstats/${meta.dbname}.pepstats"), emit: pepstats_output
    path "versions.yml", emit: versions_file

    script:
    """
    mkdir -p pepstats
    pepstats \
        -sequence ${fasta_file} \
        -outfile pepstats/${meta.dbname}.pepstats \
        -auto

    # Create versions file
    PEPSTATS_VERSION=\$(pepstats --version 2>&1 | sed 's/^EMBOSS://')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepstats: \$PEPSTATS_VERSION
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pepstats
    touch pepstats/${meta.dbname}.pepstats versions.yml
    """
}
