/*
 * BUILD_LIBRARY SUBWORKFLOW
 * Optionally run RepeatModeler to build a de novo repeat library,
 * then merge with any provided custom library.
 */

include { REPEATMODELER_REPEATMODELER } from '../modules/repeatmodeler.nf'

workflow BUILD_LIBRARY {
    take:
    genome_ch          // [ meta, fasta ]
    skip_repeatmodeler // boolean val
    custom_library     // path or null

    main:
    ch_versions = Channel.empty()
    ch_library  = Channel.empty()

    if (!skip_repeatmodeler) {
        REPEATMODELER_REPEATMODELER(genome_ch)
        ch_versions = ch_versions.mix(REPEATMODELER_REPEATMODELER.out.versions)

        // Merge RepeatModeler output with custom library if both present
        ch_library = REPEATMODELER_REPEATMODELER.out.fasta.map { meta, fa ->
            custom_library ? [meta, [fa, file(custom_library)]] : [meta, fa]
        }
    } else if (custom_library) {
        // No RepeatModeler — use custom library only
        ch_library = genome_ch.map { meta, fasta ->
            [meta, file(custom_library)]
        }
    }
    // If neither: ch_library stays empty; RepeatMasker falls back to -species

    emit:
    library  = ch_library   // [ meta, library_fa ] — may be empty
    versions = ch_versions
}
