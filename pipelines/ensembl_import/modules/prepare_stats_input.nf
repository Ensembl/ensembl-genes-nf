process PREPARE_STATS_INPUT {
    label 'process_light'
    tag "${meta.id}:statistics-input"

    input:
    tuple val(meta), path(genome_file)

    output:
    path "statistics_input_${meta.id}.csv", emit: csv
    path 'versions.yml', emit: versions

    script:
    if (!meta.db_name)
        throw new IllegalArgumentException("Sample metadata is missing db_name for ${meta.id}")
    def publishedGenome = file("${params.outdir}/refseq/${meta.id}/${genome_file.name}").toAbsolutePath().toString()

    """
    printf '%s\\n' 'dbname,species_id,genome_file' > statistics_input_${meta.id}.csv
    printf '%s\\n' "${meta.db_name},1,${publishedGenome}" >> statistics_input_${meta.id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        statistics_input: generated
    END_VERSIONS
    """

    stub:
    """
    printf '%s\\n' 'dbname,species_id,genome_file' > statistics_input_${meta.id}.csv
    printf '%s\\n' 'stub_core,1,/tmp/stub_genome.fna' >> statistics_input_${meta.id}.csv
    touch versions.yml
    """
}
