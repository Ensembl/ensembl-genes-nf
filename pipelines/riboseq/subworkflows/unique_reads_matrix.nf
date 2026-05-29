/*
 * UNIQUE READS MATRIX SUBWORKFLOW
 *
 * Builds a global count matrix from collapsed FASTA files:
 * 1. Convert collapsed FASTA -> sorted TSV (per sample)
 * 2. Build study matrices (per study, groups of samples)
 * 3. Merge into global matrix/count store (sparse Parquet by default)
 * 4. Optionally align unique reads once (single BAM for all)
 *
 * Designed for scale: 5k+ samples across 100-250 studies
 */

include { COLLAPSED_TO_TSV; COLLAPSED_TO_TSV_PARTITIONED } from '../modules/collapsed_to_tsv.nf'
include { QC_PARTITIONED_TSV } from '../modules/qc_partitioned_tsv.nf'
include { BUILD_STUDY_MATRIX; BUILD_STUDY_MATRIX_PARTITIONED } from '../modules/build_study_matrix.nf'
include { MERGE_GLOBAL_MATRIX; MERGE_GLOBAL_MATRIX_PARTITIONED } from '../modules/merge_global_matrix.nf'
include { STAR_ALIGN_UNIQUE_READS; STAR_ALIGN_UNIQUE_READS_PARTITIONED } from '../modules/star_align_unique_reads.nf'

def prefixPartitions(prefix_length) {
    def partitions = ['']
    int n = prefix_length as int
    for (int i = 0; i < n; i++) {
        partitions = partitions.collectMany { prefix -> ['A', 'C', 'G', 'T'].collect { base -> "${prefix}${base}" } }
    }
    return partitions.sort() + [('N' * n)]
}

workflow UNIQUE_READS_MATRIX {
    take:
    samples           // tuple: [ meta, collapsed_fasta ] - meta must contain 'id' and 'study_id'
    star_index        // path: STAR genome index

    main:

    if (params.matrix_use_partitioning) {
        partitions = prefixPartitions(params.matrix_partition_prefix_length ?: 4)
        partition_ordinals = Channel.fromList(
            partitions.withIndex().collect { partition, ordinal ->
                tuple(partition, ordinal, params.matrix_partition_stride ?: 1000000000L)
            }
        )

        COLLAPSED_TO_TSV_PARTITIONED(samples)
        QC_PARTITIONED_TSV(COLLAPSED_TO_TSV_PARTITIONED.out.stats)

        partitioned_tsvs = COLLAPSED_TO_TSV_PARTITIONED.out.tsvs
            .join(QC_PARTITIONED_TSV.out.qc_json)
            .flatMap { meta, tsvs, qc_json ->
                def files = tsvs instanceof List ? tsvs : [tsvs]
                files.collect { tsv ->
                    def name = tsv.name
                    def prefix = "${meta.id}."
                    def partition = name.startsWith(prefix) && name.endsWith('.tsv')
                        ? name.substring(prefix.length(), name.length() - 4)
                        : name.replaceFirst(/.*\.([A-Z]+)\.tsv$/, '$1')
                    tuple(meta.study_id ?: 'unknown_study', partition, meta.id, tsv)
                }
            }

        study_grouped_partitioned = partitioned_tsvs
            .groupTuple(by: [0, 1])
            .map { study_id, partition, sample_ids, tsv_files ->
                tuple(study_id, partition, sample_ids, tsv_files)
            }

        BUILD_STUDY_MATRIX_PARTITIONED(study_grouped_partitioned)

        partitioned_matrix_keyed = BUILD_STUDY_MATRIX_PARTITIONED.out.matrix
            .map { partition, study_id, matrix -> tuple("${partition}\t${study_id}", partition, study_id, matrix) }
        partitioned_sequences_keyed = BUILD_STUDY_MATRIX_PARTITIONED.out.sequences
            .map { partition, study_id, sequences -> tuple("${partition}\t${study_id}", sequences) }
        partitioned_metadata_keyed = BUILD_STUDY_MATRIX_PARTITIONED.out.metadata
            .map { partition, study_id, metadata -> tuple("${partition}\t${study_id}", metadata) }

        partition_study_files = partitioned_matrix_keyed
            .join(partitioned_sequences_keyed)
            .join(partitioned_metadata_keyed)
            .map { key, partition, study_id, matrix, sequences, metadata ->
                tuple(partition, [matrix, sequences, metadata])
            }
            .groupTuple(by: 0)
            .map { partition, file_groups ->
                tuple(partition, file_groups.flatten())
            }

        merge_inputs = partition_study_files
            .join(partition_ordinals)
            .map { partition, study_files, partition_ordinal, partition_stride ->
                tuple(partition, partition_ordinal, partition_stride, study_files)
            }

        MERGE_GLOBAL_MATRIX_PARTITIONED(merge_inputs)

        if (params.matrix_align_unique_reads) {
            STAR_ALIGN_UNIQUE_READS_PARTITIONED(MERGE_GLOBAL_MATRIX_PARTITIONED.out.fasta, star_index)
        }

        sample_tsvs_ch = COLLAPSED_TO_TSV_PARTITIONED.out.tsvs
        study_matrices_ch = BUILD_STUDY_MATRIX_PARTITIONED.out.matrix
        study_sequences_ch = BUILD_STUDY_MATRIX_PARTITIONED.out.sequences
        study_metadata_ch = BUILD_STUDY_MATRIX_PARTITIONED.out.metadata
        global_matrix_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.matrix
        global_counts_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.counts
        global_reads_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.reads
        global_samples_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.samples_parquet
        global_studies_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.studies_parquet
        global_retained_counts_ch = Channel.empty()
        global_fasta_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.fasta
        global_metadata_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.metadata
        global_config_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.config
        global_matrix_manifest_ch = MERGE_GLOBAL_MATRIX_PARTITIONED.out.manifest
        unique_reads_bam_ch = params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS_PARTITIONED.out.bam : Channel.empty()
        unique_reads_bai_ch = params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS_PARTITIONED.out.bai : Channel.empty()
        unique_reads_log_ch = params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS_PARTITIONED.out.log : Channel.empty()
        versions_ch = COLLAPSED_TO_TSV_PARTITIONED.out.versions.first()
            .mix(QC_PARTITIONED_TSV.out.versions.first())
            .mix(BUILD_STUDY_MATRIX_PARTITIONED.out.versions.first())
            .mix(MERGE_GLOBAL_MATRIX_PARTITIONED.out.versions)
            .mix(params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS_PARTITIONED.out.versions : Channel.empty())
    } else {
        COLLAPSED_TO_TSV(samples)

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

        study_files = BUILD_STUDY_MATRIX.out.matrix
            .join(BUILD_STUDY_MATRIX.out.sequences)
            .join(BUILD_STUDY_MATRIX.out.metadata)
            .map { study_id, matrix, sequences, metadata ->
                [matrix, sequences, metadata]
            }
            .flatten()
            .collect()

        MERGE_GLOBAL_MATRIX(study_files)

        if (params.matrix_align_unique_reads) {
            STAR_ALIGN_UNIQUE_READS(MERGE_GLOBAL_MATRIX.out.fasta, star_index)
        }

        sample_tsvs_ch = COLLAPSED_TO_TSV.out.tsv
        study_matrices_ch = BUILD_STUDY_MATRIX.out.matrix
        study_sequences_ch = BUILD_STUDY_MATRIX.out.sequences
        study_metadata_ch = BUILD_STUDY_MATRIX.out.metadata
        global_matrix_ch = MERGE_GLOBAL_MATRIX.out.matrix
        global_counts_ch = MERGE_GLOBAL_MATRIX.out.counts
        global_reads_ch = MERGE_GLOBAL_MATRIX.out.reads
        global_samples_ch = MERGE_GLOBAL_MATRIX.out.samples_parquet
        global_studies_ch = MERGE_GLOBAL_MATRIX.out.studies_parquet
        global_retained_counts_ch = MERGE_GLOBAL_MATRIX.out.retained_counts
        global_fasta_ch = MERGE_GLOBAL_MATRIX.out.fasta
        global_metadata_ch = MERGE_GLOBAL_MATRIX.out.metadata
        global_config_ch = MERGE_GLOBAL_MATRIX.out.config
        global_matrix_manifest_ch = MERGE_GLOBAL_MATRIX.out.manifest
        unique_reads_bam_ch = params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS.out.bam : Channel.empty()
        unique_reads_bai_ch = params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS.out.bai : Channel.empty()
        unique_reads_log_ch = params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS.out.log : Channel.empty()
        versions_ch = COLLAPSED_TO_TSV.out.versions.first()
            .mix(BUILD_STUDY_MATRIX.out.versions.first())
            .mix(MERGE_GLOBAL_MATRIX.out.versions)
            .mix(params.matrix_align_unique_reads ? STAR_ALIGN_UNIQUE_READS.out.versions : Channel.empty())
    }

    emit:
    // Per-sample TSVs
    sample_tsvs = sample_tsvs_ch

    // Per-study outputs
    study_matrices = study_matrices_ch
    study_sequences = study_sequences_ch
    study_metadata = study_metadata_ch

    // Global outputs
    global_matrix = global_matrix_ch
    global_counts = global_counts_ch
    global_reads = global_reads_ch
    global_samples = global_samples_ch
    global_studies = global_studies_ch
    global_retained_counts = global_retained_counts_ch
    global_fasta = global_fasta_ch
    global_metadata = global_metadata_ch
    global_config = global_config_ch
    global_matrix_manifest = global_matrix_manifest_ch

    // Alignment outputs
    unique_reads_bam = unique_reads_bam_ch
    unique_reads_bai = unique_reads_bai_ch
    unique_reads_log = unique_reads_log_ch

    // Version tracking
    versions = versions_ch
}
