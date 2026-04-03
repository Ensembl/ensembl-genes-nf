/*
 * RUN_IGTR SUBWORKFLOW
 * Split IGTR proteins into batches, run GenBlast per batch against the
 * softmasked genome, convert + filter GFF output, then cluster overlapping
 * models to select one best model per locus.
 */

include { GENBLAST        } from '../modules/genblast.nf'
include { CONVERT_GENBLAST } from '../modules/convert_genblast.nf'
include { CLUSTER_IGTR    } from '../modules/cluster_igtr.nf'

workflow RUN_IGTR {
    take:
    protein_batches   // [ meta, protein_fasta ] — one per batch (from splitFasta)
    genome            // path: softmasked genome FASTA

    main:
    ch_versions = Channel.empty()

    // Align each protein batch against the genome
    GENBLAST(protein_batches, genome)
    ch_versions = ch_versions.mix(GENBLAST.out.versions)

    // Re-attach the protein batch so convert_genblast can read biotypes
    ch_gff_with_proteins = GENBLAST.out.gff
        .join(protein_batches.map { meta, fa -> [meta, fa] })
        .map { meta, gff, proteins -> [meta, gff, proteins] }

    CONVERT_GENBLAST(
        ch_gff_with_proteins.map { meta, gff, prot -> [meta, gff] },
        ch_gff_with_proteins.map { meta, gff, prot -> prot }
    )
    ch_versions = ch_versions.mix(CONVERT_GENBLAST.out.versions)

    // Collect all filtered GFF3s under a single genome meta, then cluster
    ch_all_gff3 = CONVERT_GENBLAST.out.gff3
        .map { meta, gff3 -> gff3 }
        .collect()
        .map { gff3_list -> [[id: 'igtr'], gff3_list] }

    CLUSTER_IGTR(ch_all_gff3)
    ch_versions = ch_versions.mix(CLUSTER_IGTR.out.versions)

    emit:
    gff3     = CLUSTER_IGTR.out.gff3  // [ meta, clustered.gff3 ]
    versions = ch_versions
}
