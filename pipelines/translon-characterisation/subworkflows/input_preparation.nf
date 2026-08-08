include { NORMALISE_INTERVALS } from '../modules/normalise_intervals.nf'

workflow INPUT_PREPARATION {
    take:
    intervals
    gencode_gff3
    verdicts
    proteome
    genome

    main:
    interval_input = Channel.of([[id: 'translon-input'], intervals])
    NORMALISE_INTERVALS(interval_input)
    annotation_input = Channel.of([[id: 'gencode-reference'], gencode_gff3])
    verdict_input = Channel.of([[id: 'translonscorer-verdicts'], verdicts])
    proteome_input = Channel.of([[id: 'gencode-proteome'], proteome])
    genome_input = Channel.of([[id: 'genome-reference'], genome])

    emit:
    intervals = NORMALISE_INTERVALS.out.intervals
    annotation = annotation_input
    verdicts = verdict_input
    proteome = proteome_input
    genome = genome_input
    versions = NORMALISE_INTERVALS.out.versions
}
