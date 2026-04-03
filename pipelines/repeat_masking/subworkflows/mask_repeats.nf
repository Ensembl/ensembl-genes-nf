/*
 * MASK_REPEATS SUBWORKFLOW
 * Run RepeatMasker on genome chunks, plus RED, TRF, DUST.
 * Collect all BED outputs, merge into a single GFF3, and produce
 * a softmasked genome FASTA.
 */

include { REPEATMASKER_REPEATMASKER } from '../modules/repeatmasker.nf'
include { RED                       } from '../modules/red.nf'
include { TRF                       } from '../modules/trf.nf'
include { DUSTMASKER                } from '../modules/dustmasker.nf'
include { MERGE_REPEATS             } from '../modules/merge_repeats.nf'
include { BEDTOOLS_MASKFASTA        } from '../modules/bedtools_maskfasta.nf'

workflow MASK_REPEATS {
    take:
    chunks_ch      // [ meta, fasta ] — chunked genome slices
    genome_ch      // [ meta, fasta ] — whole genome for final masking
    library_ch     // [ meta, fasta ] — repeat library (may be empty channel)
    skip_red       // boolean
    skip_trf       // boolean
    skip_dust      // boolean

    main:
    ch_versions = Channel.empty()
    ch_beds     = Channel.empty()

    // -- RepeatMasker (run on each chunk) ------------------------------------
    // If library_ch is empty, pass [] so the module falls back to -species
    ch_lib_or_empty = library_ch.ifEmpty([[id: 'none'], []])

    ch_rm_input = chunks_ch.combine(ch_lib_or_empty.map { meta, lib -> lib })

    REPEATMASKER_REPEATMASKER(
        ch_rm_input.map { meta, fasta, lib -> [meta, fasta] },
        ch_rm_input.map { meta, fasta, lib -> lib }
    )
    ch_versions = ch_versions.mix(REPEATMASKER_REPEATMASKER.out.versions)

    // RepeatMasker .gff output (GFF2 format) is used as our BED proxy
    ch_beds = ch_beds.mix(
        REPEATMASKER_REPEATMASKER.out.gff.map { meta, gff -> [meta, gff] }
    )

    // -- RED -----------------------------------------------------------------
    if (!skip_red) {
        RED(chunks_ch)
        ch_versions = ch_versions.mix(RED.out.versions)
        ch_beds     = ch_beds.mix(RED.out.bed)
    }

    // -- TRF -----------------------------------------------------------------
    if (!skip_trf) {
        TRF(chunks_ch)
        ch_versions = ch_versions.mix(TRF.out.versions)
        ch_beds     = ch_beds.mix(TRF.out.bed)
    }

    // -- DUST ----------------------------------------------------------------
    if (!skip_dust) {
        DUSTMASKER(chunks_ch)
        ch_versions = ch_versions.mix(DUSTMASKER.out.versions)
        ch_beds     = ch_beds.mix(DUSTMASKER.out.bed)
    }

    // -- Merge all BEDs into one GFF3 ----------------------------------------
    ch_all_beds = ch_beds
        .map { meta, bed -> [meta.genome_id ?: meta.id, meta, bed] }
        .groupTuple(by: 0)
        .map { genome_id, metas, beds ->
            [[id: genome_id], beds.flatten()]
        }

    MERGE_REPEATS(ch_all_beds)
    ch_versions = ch_versions.mix(MERGE_REPEATS.out.versions)

    // -- Softmask the full genome FASTA with the merged BED ------------------
    ch_mask_input = genome_ch.join(MERGE_REPEATS.out.bed)

    BEDTOOLS_MASKFASTA(
        ch_mask_input.map { meta, fasta, bed -> [meta, bed]   },
        ch_mask_input.map { meta, fasta, bed -> [meta, fasta] }
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)

    emit:
    softmasked_fasta = BEDTOOLS_MASKFASTA.out.fasta  // [ meta, softmasked.fa ]
    repeat_gff3      = MERGE_REPEATS.out.gff3        // [ meta, repeats.gff3  ]
    repeat_bed       = MERGE_REPEATS.out.bed         // [ meta, repeats.bed   ]
    versions         = ch_versions
}
