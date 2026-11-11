nextflow.enable.dsl = 2

include { PREPARE_REFERENCES } from '../modules/prepare_references.nf'
include { BUILD_STAR_INDEX } from '../modules/build_star_index.nf'

workflow PREPARE_REFERENCES_WORKFLOW {
    take:
    references

    main:
    ch_reference = PREPARE_REFERENCES(references).reference
    ch_index = BUILD_STAR_INDEX(ch_reference).index

    emit:
    reference = ch_reference
    index = ch_index
    annotation = ch_reference.map { meta, _fasta, gtf, _chrom_sizes -> tuple(meta.assembly, gtf) }
}
