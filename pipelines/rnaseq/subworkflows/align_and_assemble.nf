/*
 * ALIGN_AND_ASSEMBLE SUBWORKFLOW
 * STAR alignment + StringTie2 assembly per sample.
 * Emits per-sample filtered GFF3.
 */

include { STAR_ALIGN      } from '../modules/star_align.nf'
include { STRINGTIE       } from '../modules/stringtie.nf'
include { FILTER_STRINGTIE } from '../modules/filter_stringtie.nf'

workflow ALIGN_AND_ASSEMBLE {
    take:
    samples    // [ meta, [ reads ] ] — meta includes .id and .strandedness
    star_index // path: STAR index directory

    main:
    ch_versions = Channel.empty()

    STAR_ALIGN(samples, star_index)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    STRINGTIE(STAR_ALIGN.out.bam)
    ch_versions = ch_versions.mix(STRINGTIE.out.versions)

    FILTER_STRINGTIE(STRINGTIE.out.gtf)
    ch_versions = ch_versions.mix(FILTER_STRINGTIE.out.versions)

    emit:
    gff3     = FILTER_STRINGTIE.out.gff3  // [ meta, sample.rnaseq.gff3 ]
    bam      = STAR_ALIGN.out.bam         // [ meta, bam ]
    sj       = STAR_ALIGN.out.sj          // [ meta, SJ.out.tab ]
    versions = ch_versions
}
