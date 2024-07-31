
include { RAPID } from '../../modules/local/rapid.nf'
include { LIST_MIRMACHINE_CLADES; FORMAT_CLADES; MATCH_CLADE } from '../../modules/local/clade_matcher.nf'
include { MIRMACHINE } from '../../modules/local/mirMachine.nf'

workflow mirMachine {

    take:
        species
        accession

    main:
        fasta_ch                =   RAPID(species, accession)
        clades_ch               =   LIST_MIRMACHINE_CLADES()
        formatted_clades_ch     =   FORMAT_CLADES(clades_ch)
        match_ch                =   MATCH_CLADE(species, accession, formatted_clades_ch)
        //match_ch.view()

        clean_species = species.replaceAll(~/\s/,"")
        //output contains several lines of the script processing, to remove all unwanted info:
        clean_nodes = match_ch.map { it ->
            it.toString().tokenize().last()
        }
        //That will give us something on the lines of: ClosestClade;Deuterostomia/Protostomia, so we need to split it:
        split_nodes = clean_nodes.map { it.split(';') }
        //and assing to each
        closer_clade = split_nodes.map { it[0].trim() }
        model = split_nodes.map { it[1].trim() }

        MIRMACHINE(clean_species, closer_clade, model, fasta_ch)

    emit:
        fasta_ch
}
