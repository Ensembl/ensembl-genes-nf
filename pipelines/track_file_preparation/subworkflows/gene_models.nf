include { VALIDATE_GTF } from '../modules/gtf_validation.nf'
include { GTF_TO_GENEPRED } from '../modules/gtf_to_genepred.nf'
include { GENEPRED_TO_BED12 } from '../modules/genepred_to_bed12.nf'
include { SORT_GENE_MODELS } from '../modules/sort_bed12.nf'
include { BED12_TO_BIGBED } from '../modules/bed12_to_bigbed.nf'
include { VALIDATE_GENE_BIGBED } from '../modules/validate_gene_bigbed.nf'
include { CHECKSUM_TRACK } from '../modules/checksum_track.nf'
include { WRITE_TRACK_RECORD } from '../modules/write_track_record.nf'

workflow PREPARE_GENE_MODEL_TRACKS {
    take:
    gene_entities // tuple [meta, gtf, chrom_sizes]

    main:
    checked = VALIDATE_GTF(gene_entities)
    genepred = GTF_TO_GENEPRED(checked.validated)
    bed12 = GENEPRED_TO_BED12(genepred.genepred)
    sorted = SORT_GENE_MODELS(bed12.bed12)
    bigbed = BED12_TO_BIGBED(sorted.sorted)
    validated = VALIDATE_GENE_BIGBED(bigbed.bigbed)
    record_input = validated.validated.map { meta, track, source -> tuple(meta + [current_track: 'gene_model'], track, source) }
    checksummed = CHECKSUM_TRACK(record_input)
    records = WRITE_TRACK_RECORD(checksummed.checksummed)

    emit:
    results = records.results
    versions = records.versions
}
