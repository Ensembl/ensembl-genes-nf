nextflow.enable.dsl = 2

include { DOWNLOAD_REFERENCES } from '../modules/download_references.nf'
include { PREPARE_REFERENCE_BUNDLE } from '../modules/prepare_reference_bundle.nf'
include { BUILD_STAR_INDEX } from '../modules/build_star_index.nf'

workflow PREPARE_REFERENCES_WORKFLOW {
    take:
    references

    main:
    ch_downloaded = DOWNLOAD_REFERENCES(references, file("${projectDir}/bin/cache_reference.py")).reference
    ch_reference = PREPARE_REFERENCE_BUNDLE(ch_downloaded, file("${projectDir}/bin/prepare_reference.py")).reference
    ch_index = BUILD_STAR_INDEX(ch_reference).index

    emit:
    reference = ch_reference
    index = ch_index
    annotation = ch_reference.map { meta, _fasta, gtf, _chrom_sizes, _rrna_bed, _excluded -> tuple(meta.assembly, gtf) }
    excluded_contigs = ch_index.map { meta, index, excluded -> tuple(meta.assembly, excluded) }
}
