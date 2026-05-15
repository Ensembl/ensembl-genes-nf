include { FETCH } from '../../modules/local/fetch'
include { MIRMACHINE } from '../../modules/local/mirmachine'

workflow mirMachine {
    take:
        input_ch
        previously_run
        fasta_dir

    main:
        genomes_to_process = input_ch.map { meta, species, accession ->
            def fasta_path = file("${fasta_dir}/${meta.id}.fa")
            if (!fasta_path.exists()) {
                return tuple(meta, species, accession)
            }
            return null
        }.filter { row -> row != null }

        rapid_ch = FETCH(genomes_to_process)

        fasta_ch = input_ch.map { meta, _species, _accession ->
            def fasta_path = file("${fasta_dir}/${meta.id}.fa")
            if (fasta_path.exists()) {
                return tuple(meta, fasta_path)
            }
            return null
        }.filter { row -> row != null }
        .mix(rapid_ch)

        mirmachine_to_run = fasta_ch
            .map { meta, fasta ->
                tuple(meta.id, meta, fasta)
            }
            .join(previously_run, by: [0], remainder:true)
            .filter { row -> row[1] != null }
            .map { _id, meta, fasta, _prev_run ->
                tuple(meta, fasta)
            }

        MIRMACHINE(mirmachine_to_run)

    emit:
        fasta = fasta_ch
        results = MIRMACHINE.out.predictions
        logs = MIRMACHINE.out.log
}
