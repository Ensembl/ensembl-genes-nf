/*
 * RUN_HOMOLOGY SUBWORKFLOW
 * Split UniProt proteins into batches, run GenBlast per batch, classify into
 * genblast_1..7 tiers, then remove redundant lower-priority overlapping models.
 */

include { GENBLAST          } from '../modules/genblast.nf'
include { CLASSIFY_GENBLAST } from '../modules/classify_genblast.nf'
include { REMOVE_REDUNDANT  } from '../modules/remove_redundant.nf'

workflow RUN_HOMOLOGY {
    take:
    protein_batches   // [ meta, protein_fasta ] — one per batch
    genome            // path: softmasked genome FASTA

    main:
    ch_versions = Channel.empty()

    GENBLAST(protein_batches, genome)
    ch_versions = ch_versions.mix(GENBLAST.out.versions)

    CLASSIFY_GENBLAST(GENBLAST.out.gff)
    ch_versions = ch_versions.mix(CLASSIFY_GENBLAST.out.versions)

    // Collect all classified GFF3s under single genome meta, then deduplicate
    ch_all_gff3 = CLASSIFY_GENBLAST.out.gff3
        .map { meta, gff3 -> gff3 }
        .collect()
        .map { gff3_list -> [[id: 'genblast_homology'], gff3_list] }

    REMOVE_REDUNDANT(ch_all_gff3)
    ch_versions = ch_versions.mix(REMOVE_REDUNDANT.out.versions)

    emit:
    gff3     = REMOVE_REDUNDANT.out.gff3  // [ meta, final.gff3 ]
    versions = ch_versions
}
