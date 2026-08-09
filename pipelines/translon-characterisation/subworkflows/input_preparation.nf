include { NORMALISE_INTERVALS } from '../modules/normalise_intervals.nf'

workflow INPUT_PREPARATION {
    take:
    intervals
    gencode_gff3
    verdicts
    proteome
    genome

    main:
    interval_input = intervals.map { meta, path -> tuple(meta, path) }
    NORMALISE_INTERVALS(interval_input)
    annotation_input = gencode_gff3.map { path -> tuple([id: 'gencode-reference'], path) }
    verdict_input = verdicts.map { meta, path -> tuple(meta + [input_type: 'translonscorer-verdicts'], path) }
    proteome_input = proteome.map { path -> tuple([id: 'gencode-proteome'], path) }
    genome_input = genome.map { path -> tuple([id: 'genome-reference'], path) }

    emit:
    intervals = NORMALISE_INTERVALS.out.intervals
    annotation = annotation_input
    verdicts = verdict_input
    proteome = proteome_input
    genome = genome_input
    versions = NORMALISE_INTERVALS.out.versions
}
