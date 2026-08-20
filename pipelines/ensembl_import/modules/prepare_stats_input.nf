process PREPARE_STATS_INPUT {
    label 'process_light'
    tag "${meta.id}:statistics-input"

    input:
    tuple val(meta), path(genome_file)

    output:
    path "statistics_input_${meta.id}.csv", emit: csv
    path 'versions.yml', emit: versions

    script:
    def speciesParts = meta.species.tokenize(' ')
    if (speciesParts.size() < 2)
        throw new IllegalArgumentException("Species must contain a genus and species: ${meta.species}")
    def speciesToken = "${speciesParts[0].toLowerCase()}_${speciesParts[1].toLowerCase()}"
    def accessionToken = meta.id.toLowerCase().replace('_', '').replaceFirst(/\./, 'v')
    def dbName = "${speciesToken}_${accessionToken}_rs_core_114_1"
    def publishedGenome = file("${params.outdir}/refseq/${meta.id}/${genome_file.name}").absolutePath

    """
    printf '%s\\n' 'taxon_id,gca,dbname,species_id,genome_file' > statistics_input_${meta.id}.csv
    printf '%s\\n' 'UNKNOWN,UNKNOWN,${dbName},1,${publishedGenome}' >> statistics_input_${meta.id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        statistics_input: generated
    END_VERSIONS
    """

    stub:
    """
    printf '%s\\n' 'taxon_id,gca,dbname,species_id,genome_file' > statistics_input_${meta.id}.csv
    printf '%s\\n' 'UNKNOWN,UNKNOWN,stub_core,1,/tmp/stub_genome.fna' >> statistics_input_${meta.id}.csv
    touch versions.yml
    """
}
