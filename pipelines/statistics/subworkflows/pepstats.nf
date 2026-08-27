include { DB_METADATA } from '../modules/db_metadata.nf'
include { FETCH_PROTEINS_ALL } from '../modules/fetch_proteins.nf'
include { RUN_PEPSTATS } from '../modules/run_pepstats.nf'
include { PARSE_PEPSTATS } from '../modules/parse_pepstats.nf'

workflow PEPSTATS {
    take:
    csvFile

    main:
    data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
        .splitCsv(sep: ',', header: true)
        .map { row ->
            if (!row.get('dbname'))
                error('Pepstats requires a dbname column in the input CSV')
            [
                gca: 'UNKNOWN',
                taxon_id: 'UNKNOWN',
                dbname: row.get('dbname'),
                species_id: row.get('species_id') ? row.get('species_id') as Integer : 1,
                protein_file: row.get('protein_file')
            ]
        }

    metadata = DB_METADATA(data).metadata.map { meta, metadata_file ->
        def values = metadata_file.text.readLines().collectEntries { line ->
            def parts = line.split('=', 2)
            [(parts[0]): parts[1]]
        }
        meta + [
            taxon_id: values.taxon_id,
            gca: values.gca,
            production_name: values.production_name
        ]
    }

    proteins = FETCH_PROTEINS_ALL(metadata).fasta_file_output
    pepstats = RUN_PEPSTATS(proteins).pepstats_output
    PARSE_PEPSTATS(pepstats)

    emit:
    proteins = proteins
    versions = DB_METADATA.out.versions_file
        .mix(FETCH_PROTEINS_ALL.out.versions_file)
        .mix(RUN_PEPSTATS.out.versions_file)
        .mix(PARSE_PEPSTATS.out.versions)
}
