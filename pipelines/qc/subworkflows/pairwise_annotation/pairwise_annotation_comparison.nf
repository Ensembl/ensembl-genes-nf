include { STAGE_ANNOTATION_PAIR } from '../../modules/pairwise_annotation/stage_annotation_pair.nf'
include { PAIR_GENES } from '../../modules/pairwise_annotation/pair_genes.nf'
include { PAIRWISE_TRANSCRIPT_CONCORDANCE } from '../../modules/pairwise_annotation/transcript_concordance.nf'
include { PAIRWISE_CODING_INTEGRITY } from '../../modules/pairwise_annotation/coding_integrity.nf'
include { PAIRWISE_GENE_PRESENCE } from '../../modules/pairwise_annotation/gene_presence.nf'
include { PAIRWISE_MULTI_MAPPING } from '../../modules/pairwise_annotation/multi_mapping.nf'
include { PAIRWISE_GFF_FEATURE_METRICS } from '../../modules/pairwise_annotation/gff_feature_metrics.nf'
include { PAIRWISE_COUNT_GFF_TRANSCRIPTS } from '../../modules/pairwise_annotation/count_gff_transcripts.nf'
include { PAIRWISE_GFF_TO_GTF as SOURCE_A_GFF_TO_GTF } from '../../modules/pairwise_annotation/gff_to_gtf.nf'
include { PAIRWISE_GFF_TO_GTF as SOURCE_B_GFF_TO_GTF } from '../../modules/pairwise_annotation/gff_to_gtf.nf'
include { PAIRWISE_GFFCOMPARE; PAIRWISE_PARSE_GFFCOMPARE } from '../../modules/pairwise_annotation/gffcompare.nf'

workflow PAIRWISE_ANNOTATION_COMPARISON {

    take:
        pairwise_samples_ch
        ensg_lookup

    main:
        STAGE_ANNOTATION_PAIR(pairwise_samples_ch, params.pairwise_cache_dir)

        staged = STAGE_ANNOTATION_PAIR.out.staged

        PAIR_GENES(staged)

        staged_gffs = staged.map { meta, source_a_gff, source_b_gff, assembly_report ->
            tuple(meta, source_a_gff, source_b_gff)
        }

        rbh_inputs = staged_gffs.join(PAIR_GENES.out.rbh)
        all_pair_inputs = PAIR_GENES.out.all_pairs

        if (params.run_pairwise_transcript_concordance) {
            PAIRWISE_TRANSCRIPT_CONCORDANCE(rbh_inputs)
            transcript_concordance_ch = PAIRWISE_TRANSCRIPT_CONCORDANCE.out.metrics
        }
        else {
            transcript_concordance_ch = Channel.empty()
        }

        if (params.run_pairwise_coding_integrity) {
            PAIRWISE_CODING_INTEGRITY(rbh_inputs)
            coding_integrity_ch = PAIRWISE_CODING_INTEGRITY.out.metrics
        }
        else {
            coding_integrity_ch = Channel.empty()
        }

        if (params.run_pairwise_gene_presence) {
            PAIRWISE_GENE_PRESENCE(staged_gffs, ensg_lookup)
            gene_presence_ch = PAIRWISE_GENE_PRESENCE.out.metrics
        }
        else {
            gene_presence_ch = Channel.empty()
        }

        if (params.run_pairwise_multi_mapping) {
            PAIRWISE_MULTI_MAPPING(all_pair_inputs)
            multi_mapping_ch = PAIRWISE_MULTI_MAPPING.out.metrics
        }
        else {
            multi_mapping_ch = Channel.empty()
        }

        if (params.run_pairwise_gff_feature_metrics) {
            PAIRWISE_GFF_FEATURE_METRICS(staged_gffs)
            feature_counts_ch = PAIRWISE_GFF_FEATURE_METRICS.out.features
            gene_metrics_ch = PAIRWISE_GFF_FEATURE_METRICS.out.gene_metrics
            tx_metrics_ch = PAIRWISE_GFF_FEATURE_METRICS.out.tx_metrics
        }
        else {
            feature_counts_ch = Channel.empty()
            gene_metrics_ch = Channel.empty()
            tx_metrics_ch = Channel.empty()
        }

        if (params.run_pairwise_transcript_counts) {
            transcript_count_inputs = staged
                .flatMap { meta, source_a_gff, source_b_gff, assembly_report ->
                    [
                        tuple(meta, 'source_a', source_a_gff),
                        tuple(meta, 'source_b', source_b_gff)
                    ]
                }

            PAIRWISE_COUNT_GFF_TRANSCRIPTS(transcript_count_inputs)
            transcript_counts_ch = PAIRWISE_COUNT_GFF_TRANSCRIPTS.out.counts
        }
        else {
            transcript_counts_ch = Channel.empty()
        }

        if (params.run_pairwise_gffcompare) {
            SOURCE_A_GFF_TO_GTF(
                staged.map { meta, source_a_gff, source_b_gff, assembly_report ->
                    tuple(meta, 'source_a', source_a_gff)
                }
            )
            SOURCE_B_GFF_TO_GTF(
                staged.map { meta, source_a_gff, source_b_gff, assembly_report ->
                    tuple(meta, 'source_b', source_b_gff)
                }
            )

            SOURCE_A_GFF_TO_GTF.out.gtf
                .join(SOURCE_B_GFF_TO_GTF.out.gtf)
                .flatMap { meta, source_a_label, source_a_gtf, source_b_label, source_b_gtf ->
                    [
                        tuple(meta, 'source_a_to_source_b', source_a_gtf, source_b_gtf),
                        tuple(meta, 'source_b_to_source_a', source_b_gtf, source_a_gtf)
                    ]
                }
                .set { gffcompare_inputs }

            PAIRWISE_GFFCOMPARE(gffcompare_inputs)
            PAIRWISE_PARSE_GFFCOMPARE(PAIRWISE_GFFCOMPARE.out.tmap)

            gffcompare_tmap_ch = PAIRWISE_GFFCOMPARE.out.tmap
            gffcompare_stats_ch = PAIRWISE_GFFCOMPARE.out.stats
            gffcompare_class_counts_ch = PAIRWISE_PARSE_GFFCOMPARE.out.counts
        }
        else {
            gffcompare_tmap_ch = Channel.empty()
            gffcompare_stats_ch = Channel.empty()
            gffcompare_class_counts_ch = Channel.empty()
        }

    emit:
        staged_inputs = staged
        rbh = PAIR_GENES.out.rbh
        all_pairs = PAIR_GENES.out.all_pairs
        assembly_summary = PAIR_GENES.out.summary
        transcript_concordance = transcript_concordance_ch
        coding_integrity = coding_integrity_ch
        gene_presence = gene_presence_ch
        multi_mapping = multi_mapping_ch
        feature_counts = feature_counts_ch
        gene_metrics = gene_metrics_ch
        tx_metrics = tx_metrics_ch
        transcript_counts = transcript_counts_ch
        gffcompare_tmap = gffcompare_tmap_ch
        gffcompare_stats = gffcompare_stats_ch
        gffcompare_class_counts = gffcompare_class_counts_ch
        versions = STAGE_ANNOTATION_PAIR.out.versions
            .mix(PAIR_GENES.out.versions)
}
