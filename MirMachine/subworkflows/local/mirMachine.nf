include { RAPID } from '../../modules/local/rapid'
include { LIST_MIRMACHINE_CLADES } from '../../modules/local/list_mirmachine_clades'
include { FORMAT_CLADES } from '../../modules/local/format_clades'
include { MATCH_CLADE } from '../../modules/local/match_clade'
include { MIRMACHINE } from '../../modules/local/mirmachine'

workflow mirMachine {
    take:
        input_ch
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
        rapid_ch = RAPID(genomes_to_process)

        // Combine existing and newly downloaded FASTA files
        fasta_ch = input_ch.map { meta, species, accession ->
            def fasta_path = file("${fasta_dir}/${meta.id}.fa")
            if (fasta_path.exists()) {
                return tuple(meta, fasta_path)
            }
            return null
        }.filter { it != null }
        .mix(rapid_ch)

        clades_ch = LIST_MIRMACHINE_CLADES()
        formatted_clades_ch = FORMAT_CLADES(clades_ch.clades)

        // Match clades
        match_ch = MATCH_CLADE(fasta_ch, formatted_clades_ch.formatted_clades)

        // Prepare input for MIRMACHINE
        mirmachine_input = fasta_ch
            .join(match_ch, by: [0])
            .map { meta, fasta, clade_info_file ->
                def clade_info = clade_info_file.text.trim()
                def (closer_clade, model) = clade_info.tokenize(';')
                tuple(meta, fasta, closer_clade.trim(), model.trim())
            }

        // Run MIRMACHINE
        MIRMACHINE(mirmachine_input)

    emit:
        fasta = fasta_ch
        results = MIRMACHINE.out.predictions
        logs = MIRMACHINE.out.log
}