include { FETCH } from '../../modules/local/fetch'
include { GET_NCBI_TAXONOMY } from '../../modules/local/get_ncbi_taxonomy'
include { LIST_MIRMACHINE_CLADES } from '../../modules/local/list_mirmachine_clades'
include { FORMAT_CLADES } from '../../modules/local/format_clades'
include { MATCH_CLADE } from '../../modules/local/match_clade'
include { MIRMACHINE } from '../../modules/local/mirmachine'
include { COLLATE_RESULTS } from '../../modules/local/collate_results'
include { GENERATE_SCORES } from '../../modules/local/generate_scores'

workflow mirMachine {
    take:
        input_ch
        previously_run
        fasta_dir

    main:
        // Pre-flight check: determine which genomes need to be processed
        genomes_to_process = input_ch.map { meta, species, accession ->
            def fasta_path = file("${fasta_dir}/${meta.id}.fa")
            if (!fasta_path.exists()) {
                return tuple(meta, species, accession)
            }
            return null
        }.filter { it != null }

        // Download missing FASTA files
        rapid_ch = FETCH(genomes_to_process)

        // Combine existing and newly downloaded FASTA files
        fasta_ch = input_ch.map { meta, species, accession ->
            def fasta_path = file("${fasta_dir}/${meta.id}.fa")
            if (fasta_path.exists()) {
                return tuple(meta, fasta_path)
            }
            return null
        }.filter { it != null }
        .mix(rapid_ch)

        taxa_file = file("$params.outdir/taxa/taxa.sqlite")

        if (taxa_file.exists()) {
            taxa_ch = Channel.fromPath("$params.outdir/taxa/taxa.sqlite").first()
        } else {
            taxa_ch = GET_NCBI_TAXONOMY()
        }

        clades_ch = LIST_MIRMACHINE_CLADES()
        formatted_clades_ch = FORMAT_CLADES(clades_ch.clades)
        
	    // Match clades
        match_ch = MATCH_CLADE(fasta_ch, formatted_clades_ch.formatted_clades, taxa_ch)

        // Prepare input for MIRMACHINE
        mirmachine_input = fasta_ch
            .join(match_ch, by: [0])
            .map { meta, fasta, clade_info_file ->
                def clade_info = clade_info_file.text.trim()
                def (closer_clade, model) = clade_info.tokenize(';')
                tuple(meta.id, meta, fasta, closer_clade.trim(), model.trim())
            }

        mirmachine_to_run = mirmachine_input
            .join(previously_run, by: [0], remainder:true)
            .filter { it[1] != null }
            .map { id, meta, fasta, closer_clade, model, prev_run ->
                tuple(meta, fasta, closer_clade, model)
            }

        MIRMACHINE(mirmachine_to_run)

    emit:
        fasta = fasta_ch
        results = MIRMACHINE.out.predictions
        logs = MIRMACHINE.out.log
}
