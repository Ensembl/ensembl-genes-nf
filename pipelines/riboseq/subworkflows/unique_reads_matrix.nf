/*
 * UNIQUE READS MATRIX SUBWORKFLOW
 *
 * Builds a global count matrix from collapsed FASTA files:
 * 1. Convert collapsed FASTA -> sorted TSV (per sample)
 * 2. Build study matrices (per study, groups of samples)
 * 3. Merge into global matrix (Zarr format)
 * 4. Align unique reads once (single BAM for all)
 *
 * Designed for scale: 5k+ samples across 100-250 studies
 */

include { COLLAPSED_TO_TSV } from '../modules/collapsed_to_tsv.nf'
include { BUILD_STUDY_MATRIX } from '../modules/build_study_matrix.nf'
include { MERGE_GLOBAL_MATRIX } from '../modules/merge_global_matrix.nf'
include { STAR_ALIGN_UNIQUE_READS } from '../modules/star_align_unique_reads.nf'

workflow UNIQUE_READS_MATRIX {
    take:
    samples           // tuple: [ meta, collapsed_fasta ] - meta must contain 'id' and 'study_id'
    star_index        // path: STAR genome index

    main:

    //
    // STAGE 1: Convert collapsed FASTA to sorted TSV (per sample)
    //
    COLLAPSED_TO_TSV(samples)

    //
    // STAGE 2: Group by study and build study matrices
    //
    // Group TSVs by study_id, collecting sample IDs and TSV files
    study_grouped = COLLAPSED_TO_TSV.out.tsv
        .map { meta, tsv ->
            def study_id = meta.study_id ?: 'unknown_study'
            tuple(study_id, meta.id, tsv)
        }
        .groupTuple(by: 0)
        .map { study_id, sample_ids, tsv_files ->
            tuple(study_id, sample_ids, tsv_files)
        }

    study_grouped.view { study_id, sample_ids, tsv_files ->
        println "Study: ${study_id}, Samples: ${sample_ids.size()}, TSVs: ${tsv_files.size()}"
    }    

    BUILD_STUDY_MATRIX(study_grouped)

    //
    // STAGE 3: Collect all study outputs and merge into global matrix
    //
    // Collect all study files (staged flat, script reorganizes by study_id prefix)
    // Files named: {study_id}_matrix.npz, {study_id}_vocab.pkl, etc.
    study_files = BUILD_STUDY_MATRIX.out.matrix
        .join(BUILD_STUDY_MATRIX.out.vocab)
        .join(BUILD_STUDY_MATRIX.out.sequences)
        .join(BUILD_STUDY_MATRIX.out.metadata)
        .map { study_id, matrix, vocab, sequences, metadata ->
            [matrix, vocab, sequences, metadata]
        }
        .flatten()
        .collect()

    MERGE_GLOBAL_MATRIX(study_files)

    //
    // STAGE 4: Align unique reads FASTA
    //
    STAR_ALIGN_UNIQUE_READS(
        MERGE_GLOBAL_MATRIX.out.fasta,
        star_index
    )

    emit:
    // Per-sample TSVs
    sample_tsvs = COLLAPSED_TO_TSV.out.tsv              // tuple: [ meta, tsv ]

    // Per-study outputs
    study_matrices = BUILD_STUDY_MATRIX.out.matrix      // tuple: [ study_id, matrix.npz ]
    study_vocab = BUILD_STUDY_MATRIX.out.vocab          // tuple: [ study_id, vocab.pkl ]
    study_sequences = BUILD_STUDY_MATRIX.out.sequences  // tuple: [ study_id, sequences.txt.gz ]
    study_metadata = BUILD_STUDY_MATRIX.out.metadata    // tuple: [ study_id, metadata.json ]

    // Global outputs
    global_matrix = MERGE_GLOBAL_MATRIX.out.matrix      // path: global_matrix.zarr
    global_fasta = MERGE_GLOBAL_MATRIX.out.fasta        // path: unique_reads.fasta
    global_metadata = MERGE_GLOBAL_MATRIX.out.metadata  // path: read_metadata.parquet
    global_config = MERGE_GLOBAL_MATRIX.out.config      // path: index_config.json

    // Alignment outputs
    unique_reads_bam = STAR_ALIGN_UNIQUE_READS.out.bam  // path: unique_reads.bam
    unique_reads_bai = STAR_ALIGN_UNIQUE_READS.out.bai  // path: unique_reads.bam.bai

    // Version tracking
    versions = COLLAPSED_TO_TSV.out.versions.first()
        .mix(BUILD_STUDY_MATRIX.out.versions.first())
        .mix(MERGE_GLOBAL_MATRIX.out.versions)
        .mix(STAR_ALIGN_UNIQUE_READS.out.versions)
}
