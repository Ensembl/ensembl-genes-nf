include { LOAD_RIBOSEQ_OUTPUTS } from '../subworkflows/input_contract.nf'
include { MERGE_RIBO_INPUTS } from '../modules/merge_bams.nf'
include { INFLATE_UNIQUE_BAM } from '../modules/inputs/inflate_unique_bam.nf'
include { PREPARE_CALLER_INPUTS } from '../subworkflows/caller_inputs.nf'
include { PREPARE_RIBOTRICER_OFFSETS } from '../modules/preparation/ribotricer_offsets.nf'
include { PREPARE_RIBOTRICER_AUTO_OFFSETS } from '../modules/preparation/ribotricer_auto_offsets.nf'
include { MAKE_TRANSCRIPTOME_ANNOTATION } from '../modules/preparation/transcriptome_annotation.nf'
include { RUN_RIBOCODE } from '../modules/callers/ribocode.nf'
include { PREPARE_RIBOTRICER_ORFS } from '../modules/callers/ribotricer_prepare.nf'
include { RUN_RIBOTRICER } from '../modules/callers/ribotricer_detect.nf'
include { RUN_ORFQUANT } from '../modules/callers/orfquant.nf'
include { PREPARE_RPBP_GENOME } from '../modules/callers/rpbp_prepare.nf'
include { RUN_RPBP } from '../modules/callers/rpbp_run.nf'
include { IRIBO_GET_CANDIDATES } from '../modules/callers/iribo_candidates.nf'
include { IRIBO_GENERATE_PROFILE } from '../modules/callers/iribo_profile.nf'
include { IRIBO_GENERATE_TRANSLATOME } from '../modules/callers/iribo_translatome.nf'
include { MAKE_ORFRATER_TFAMS } from '../modules/callers/orfrater_make_tfams.nf'
include { FIND_ORFRATER_ORFS } from '../modules/callers/orfrater_find_orfs.nf'
include { REGRESS_ORFRATER } from '../modules/callers/orfrater_regress.nf'
include { RATE_ORFRATER } from '../modules/callers/orfrater_rate.nf'
include { RUN_ORFRATER } from '../modules/callers/orfrater_quantify.nf'
include { RUN_RIBORF } from '../modules/callers/riborf.nf'
include { RUN_RIBOTISH } from '../modules/callers/ribotish.nf'
include { PREPARE_RIBOTIE_DATA } from '../modules/callers/ribotie_prepare.nf'
include { RUN_RIBOTIE } from '../modules/callers/ribotie.nf'
include { STANDARDISE_CALLER as STANDARDISE_RIBOCODE } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_RIBOTRICER } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_ORFQUANT } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_RPBP } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_IRIBO } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_ORFRATER } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_PRICE } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_RIBORF } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_RIBOTISH } from '../modules/standardise/standardise_caller.nf'
include { STANDARDISE_CALLER as STANDARDISE_RIBOTIE } from '../modules/standardise/standardise_caller.nf'
include { GEDI_INDEXGENOME } from '../../../modules/nf-core/gedi/indexgenome/main.nf'
include { GEDI_PRICE } from '../../../modules/nf-core/gedi/price/main.nf'
include { COLLECT_ORF_CALLS } from '../modules/collect_orf_calls.nf'
include { TRANSLON_CONSENSUS_BRIDGE } from './translon_consensus_bridge.nf'
include { TRANSLON_CHARACTERISATION } from '../../translon-characterisation/workflows/translon_characterisation.nf'
include { RIBOMETRIC } from '../../riboseq/modules/ribometric.nf'
include { MAKE_PARTITION_MANIFEST } from '../modules/preparation/partition_manifest.nf'
include { PREPARE_RIBOCODE_SHARD } from '../modules/preparation/ribocode_shard.nf'
include { PREPARE_IRIBO_SHARD } from '../modules/preparation/prepare_iribo_shard.nf'
include { MERGE_IRIBO_SHARDS } from '../modules/preparation/merge_iribo_shards.nf'
include { MAKE_PRICE_CONTIG_MANIFEST; PREPARE_PRICE_CONTIG } from '../modules/preparation/price_contig.nf'
include { GEDI_PRICE_SHARDED } from '../modules/callers/gedi_price_sharded.nf'

def tool_selected(selected_tools, name) {
    def tool_groups = [
        // These sets intentionally overlap: they represent method
        // comparisons, not mutually exclusive implementation generations.
        periodicity:       ['ribocode', 'ribotricer', 'rpbp'],
        frame_tests:       ['ribocode', 'ribotish'],
        learned_models:    ['orfrater', 'riborf', 'ribotie'],
        probabilistic:     ['rpbp', 'price'],
        candidate_scoring: ['iribo', 'orfquant']
    ]
    selected_tools.contains(name) ||
    selected_tools.contains('all') ||
    selected_tools.any { selector -> tool_groups[selector]?.contains(name) }
}

workflow TRANSLON_ANALYSIS {
    ribotie_gpu_enabled = params.ribotie_gpu instanceof Boolean
        ? params.ribotie_gpu
        : params.ribotie_gpu.toString().toBoolean()

    LOAD_RIBOSEQ_OUTPUTS(
        params.riboseq_outdir,
        params.samplesheet,
        params.transcriptome_bam_glob,
        params.genome_bam_glob,
        params.offsets_glob,
        params.translonscorer_glob,
        params.merge_group
    )

    gtf = channel.value(file(params.gtf, checkIfExists: true))
    orf_gtf = params.canonical_gtf
        ? channel.value(file(params.canonical_gtf, checkIfExists: true))
        : gtf
    fasta = channel.value(file(params.fasta, checkIfExists: true))
    ribosomal_fasta = params.ribosomal_fasta ? channel.value(file(params.ribosomal_fasta, checkIfExists: true)) : channel.empty()
    adapter_fasta = params.adapter_fasta ? channel.value(file(params.adapter_fasta, checkIfExists: true)) : channel.empty()
    MERGE_RIBO_INPUTS(
        LOAD_RIBOSEQ_OUTPUTS.out.transcriptome,
        LOAD_RIBOSEQ_OUTPUTS.out.genome,
        LOAD_RIBOSEQ_OUTPUTS.out.offsets,
        params.merge_inputs,
        params.merge_group
    )
    // Input BAMs contain one record per unique sequence, with abundance in
    // the read-name suffix (_xN). Inflate both coordinate systems before any
    // caller or pooled QC sees the alignments.
    inflated_inputs = MERGE_RIBO_INPUTS.out.transcriptome
        .mix(MERGE_RIBO_INPUTS.out.genome)
    INFLATE_UNIQUE_BAM(inflated_inputs)
    tx = INFLATE_UNIQUE_BAM.out.inflated
        .filter { meta, _bam, _bai -> meta.bam_type == 'transcriptome' }
    gn = INFLATE_UNIQUE_BAM.out.inflated
        .filter { meta, _bam, _bai -> meta.bam_type == 'genome' }
    ribo_fastq = LOAD_RIBOSEQ_OUTPUTS.out.ribo_fastq
    if (params.merge_inputs) {
        if (!params.ribometric_annotation) {
            error 'Merged inputs require --ribometric_annotation so offsets can be recalculated on the pooled BAM'
        }
        ribometric_annotation = channel.value(file(params.ribometric_annotation, checkIfExists: true))
        // A pooled BAM needs pooled QC and offsets. Do not combine offsets from
        // the individual libraries: those offsets are library-specific.
        RIBOMETRIC(
            tx,
            ribometric_annotation,
            channel.value(file('NO_OFFSET_FILE'))
        )
        offsets = RIBOMETRIC.out.offsets
    } else {
        // In per-sample mode retain the QC-selected offsets produced upstream.
        offsets = MERGE_RIBO_INPUTS.out.offsets
    }
    selected_tools = params.tools.split(',').collect { tool -> tool.trim().toLowerCase() }.findAll { tool -> tool }
    start_codons = (params.start_codons ?: 'ATG').split(',').collect { codon -> codon.trim().toUpperCase() }.findAll { codon -> codon ==~ /[ACGT]{3}/ }
    if (!start_codons) error 'start_codons must contain at least one three-base DNA codon'
    run_riborf = tool_selected(selected_tools, 'riborf')
    run_orfrater = tool_selected(selected_tools, 'orfrater')
    run_rpbp = tool_selected(selected_tools, 'rpbp') && !params.skip_fastq_tools
    run_ribotie = tool_selected(selected_tools, 'ribotie')
    run_iribo = false
    ribocode_enabled = false
    // The upstream transcriptome BAM is already the native input for
    // Ribotricer, RiboTIE and most learned callers. Prepare legacy models
    // only for callers that require genePred/BED/SAM representations.
    transcript_models = channel.empty()
    riborf_inputs = channel.empty()
    ribotricer_inputs = channel.empty()
    if (tool_selected(selected_tools, 'ribotricer')) {
        if (params.ribotricer_external_offsets) {
            PREPARE_RIBOTRICER_OFFSETS(offsets)
            ribotricer_inputs = tx
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(PREPARE_RIBOTRICER_OFFSETS.out.prepared.map { meta, read_lengths, psite_offsets -> tuple(meta.id, read_lengths, psite_offsets) }, by: 0)
                .map { _id, meta, bam, bai, read_lengths, psite_offsets -> tuple(meta, bam, bai, read_lengths, psite_offsets) }
        } else {
            PREPARE_RIBOTRICER_AUTO_OFFSETS(tx)
            ribotricer_inputs = tx
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(PREPARE_RIBOTRICER_AUTO_OFFSETS.out.prepared.map { meta, read_lengths, psite_offsets -> tuple(meta.id, read_lengths, psite_offsets) }, by: 0)
                .map { _id, meta, bam, bai, read_lengths, psite_offsets -> tuple(meta, bam, bai, read_lengths, psite_offsets) }
        }
    }
    if (run_orfrater || run_riborf) {
        if (run_riborf && !params.samplesheet) {
            log.warn 'RibORF was selected but no samplesheet was supplied; skipping RibORF while continuing other callers'
            run_riborf = false
        }
        // Do not call ifEmpty here. In merged mode this channel is produced by
        // RiboMetric, so checking it during workflow construction races the
        // upstream process and can falsely report that offsets are missing.
        caller_offsets = run_riborf ? offsets : channel.empty()
        PREPARE_CALLER_INPUTS(tx, orf_gtf, caller_offsets)
        transcript_models = PREPARE_CALLER_INPUTS.out.transcript_models
        riborf_inputs = PREPARE_CALLER_INPUTS.out.riborf
    }
    // Published callers currently use one alignment contract per sample.
    // Prefer genome-space BAMs; fall back to transcriptome BAMs when necessary.
    published_inputs = gn.ifEmpty(tx)

    if (tool_selected(selected_tools, 'ribocode')) {
        MAKE_TRANSCRIPTOME_ANNOTATION(tx, orf_gtf, fasta)
        ribocode_annotation = MAKE_TRANSCRIPTOME_ANNOTATION.out.annotation
        if ((params.partition_count ?: 1) > 1 && (params.partition_mode ?: 'transcriptome') == 'transcriptome') {
            // RiboCode operates in transcriptome space, so shard it only when
            // transcriptome partitioning is requested. Genome sharding is for
            // callers such as iRibo and PRICE; RiboCode remains unsharded there.
            ribocode_enabled = true
            MAKE_PARTITION_MANIFEST(ribocode_annotation.map { meta, bam, _bai, tx_gtf, _tx_fasta -> tuple(meta, tx_gtf, bam) })
            partition_rows = MAKE_PARTITION_MANIFEST.out.manifest
                .map { meta, manifest -> tuple(meta.id, meta, manifest) }
                .flatMap { id, meta, manifest -> manifest.splitCsv(header: true, sep: '\t').collect { row -> tuple(id, meta, row) } }
            ribocode_shard_inputs = ribocode_annotation
                .map { meta, bam, bai, tx_gtf, tx_fasta -> tuple(meta.id, meta, bam, bai, tx_gtf, tx_fasta) }
                .join(partition_rows, by: 0)
                .map { _id, meta, bam, bai, tx_gtf, tx_fasta, _manifest_meta, partition ->
                    tuple(meta + [shard_id: partition.partition_id.toString()], bam, bai, tx_gtf, tx_fasta, partition)
                }
            PREPARE_RIBOCODE_SHARD(ribocode_shard_inputs)
            ribocode_inputs = PREPARE_RIBOCODE_SHARD.out.shard.flatMap { meta, bam, bai, tx_gtf, tx_fasta ->
                start_codons.collect { codon -> tuple(meta + [codon: codon], bam, bai, tx_gtf, tx_fasta, codon) }
            }
        } else {
            // Genome sharding does not apply to RiboCode. Keep this caller
            // enabled and run one transcriptome-space job per codon.
            ribocode_enabled = true
            ribocode_inputs = ribocode_annotation.flatMap { meta, bam, bai, tx_gtf, tx_fasta ->
                start_codons.collect { codon -> tuple(meta + [codon: codon, shard_id: 'all'], bam, bai, tx_gtf, tx_fasta, codon) }
            }
        }
        if (ribocode_enabled) {
            RUN_RIBOCODE(ribocode_inputs)
            STANDARDISE_RIBOCODE(RUN_RIBOCODE.out.raw, 'ribocode', orf_gtf)
        }
    }
    if (tool_selected(selected_tools, 'ribotricer')) {
        ribotricer_per_codon = ribotricer_inputs.flatMap { meta, bam, bai, read_lengths, psite_offsets ->
            start_codons.collect { codon -> tuple(meta + [codon: codon, shard_id: 'all'], bam, bai, read_lengths, psite_offsets, codon) }
        }
        PREPARE_RIBOTRICER_ORFS(ribotricer_per_codon, orf_gtf, fasta)
        RUN_RIBOTRICER(PREPARE_RIBOTRICER_ORFS.out.index)
        STANDARDISE_RIBOTRICER(RUN_RIBOTRICER.out.raw, 'ribotricer', orf_gtf)
    }
    if (tool_selected(selected_tools, 'orfquant')) {
        orfquant_alignments = gn.ifEmpty(tx)
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .join(offsets.map { meta, offset -> tuple(meta.id, offset) }, by: 0)
            .map { _id, meta, bam, bai, offset -> tuple(meta, bam, bai, offset) }
        RUN_ORFQUANT(orfquant_alignments, orf_gtf, fasta)
        STANDARDISE_ORFQUANT(RUN_ORFQUANT.out.raw, 'orfquant', orf_gtf)
    }
    if (run_rpbp) {
        if (!params.ribosomal_fasta || !params.adapter_fasta) {
            log.warn 'Rp-Bp was selected but --ribosomal_fasta and/or --adapter_fasta were not provided; skipping Rp-Bp while continuing other callers'
            run_rpbp = false
        } else {
            PREPARE_RPBP_GENOME(ribo_fastq, orf_gtf, fasta, ribosomal_fasta, adapter_fasta)
            RUN_RPBP(PREPARE_RPBP_GENOME.out.config)
            STANDARDISE_RPBP(RUN_RPBP.out.raw, 'rpbp', orf_gtf)
        }
    }
    run_iribo = tool_selected(selected_tools, 'iribo')
    if (run_iribo) {
        if ((params.partition_count ?: 1) > 1) {
            if ((params.partition_mode ?: 'genome') != 'genome' || !params.partition_fai) {
                log.warn 'iRibo sharding requires --partition_mode genome and --partition_fai; skipping iRibo while continuing other callers'
                run_iribo = false
            } else {
                run_iribo = true
            }
            if (run_iribo) {
            partition_fai = file(params.partition_fai, checkIfExists: true)
            MAKE_PARTITION_MANIFEST(published_inputs.map { values -> tuple(values[0], partition_fai, values[1]) })
            iribo_partition_manifests = MAKE_PARTITION_MANIFEST.out.manifest
                .map { meta, manifest -> tuple(meta.id, manifest) }
            iribo_shard_inputs = published_inputs
                .map { values -> tuple(values[0].id, values[0], values[1], values[2]) }
                .combine(orf_gtf)
                .combine(fasta)
                .combine(iribo_partition_manifests, by: 0)
                .flatMap { _id, meta, bam, bai, shard_gtf, shard_fasta, manifest ->
                    manifest.splitCsv(header: true, sep: '\t').collect { row ->
                        def partition = [
                            partition_id: row.partition_id as Integer,
                            mode: row.mode.toString(),
                            contig: row.contig.toString(),
                            start: row.start as Integer,
                            end: row.end as Integer,
                            padding: row.padding as Integer,
                            annotation_load: row.annotation_load as Integer,
                            estimated_read_load: row.estimated_read_load?.toString() ?: ''
                        ]
                        tuple(meta + [shard_id: partition.partition_id.toString()], bam, bai, shard_gtf, shard_fasta, partition)
                    }
                }
            PREPARE_IRIBO_SHARD(iribo_shard_inputs)
            IRIBO_GET_CANDIDATES(PREPARE_IRIBO_SHARD.out.shard.map { meta, bam, bai, shard_gtf, shard_fasta -> tuple(meta, bam, bai, shard_gtf, shard_fasta) })
            IRIBO_GENERATE_PROFILE(IRIBO_GET_CANDIDATES.out.candidates.map { meta, bam, bai, candidates, shard_gtf, shard_fasta -> tuple(meta, bam, bai, candidates, shard_gtf, shard_fasta) })
            iribo_profiles = IRIBO_GENERATE_PROFILE.out.profile
                .map { meta, profile, candidates -> tuple(meta.id, meta, profile, candidates) }
                .groupTuple(by: 0)
                .map { _id, metas, profiles, candidates -> tuple(metas[0], profiles, candidates) }
            MERGE_IRIBO_SHARDS(iribo_profiles)
            IRIBO_GENERATE_TRANSLATOME(MERGE_IRIBO_SHARDS.out.merged)
            }
        } else {
            IRIBO_GET_CANDIDATES(published_inputs
                .combine(orf_gtf)
                .combine(fasta)
                .map { meta, bam, bai, caller_gtf, caller_fasta -> tuple(meta, bam, bai, caller_gtf, caller_fasta) })
            IRIBO_GENERATE_PROFILE(IRIBO_GET_CANDIDATES.out.candidates)
            IRIBO_GENERATE_TRANSLATOME(IRIBO_GENERATE_PROFILE.out.profile)
        }
        if (run_iribo) {
            STANDARDISE_IRIBO(IRIBO_GENERATE_TRANSLATOME.out.raw, 'iribo', orf_gtf)
        }
    }
    if (run_orfrater) {
        orfrater_inputs = transcript_models
            .map { meta, bam, bai, _genepred, bed12 -> tuple(meta, bam, bai, bed12) }
        if (params.orfrater_model) {
            orfrater_model_path = file(params.orfrater_model, checkIfExists: true)
            required_orfrater_files = ['orfratings.h5', 'metagene.txt', 'offsets.txt']
            missing_orfrater_files = required_orfrater_files.findAll { name -> !file("${orfrater_model_path}/${name}").exists() }
            if (missing_orfrater_files) {
                log.warn "ORF-RATER model directory is missing: ${missing_orfrater_files.join(', ')} (${orfrater_model_path}); skipping ORF-RATER while continuing other callers"
                run_orfrater = false
            } else {
                orfrater_inputs_with_model = orfrater_inputs
                    .map { meta, bam, bai, bed12 -> tuple(meta, bam, bai, bed12, orfrater_model_path) }
                RUN_ORFRATER(orfrater_inputs_with_model)
            }
        } else {
            training_inputs = transcript_models
                .map { meta, bam, bai, _genepred, bed12 -> tuple(meta.id, meta, bam, bai, bed12) }
                .join(offsets.map { meta, offset -> tuple(meta.id, offset) }, by: 0)
                .map { _id, meta, bam, bai, bed12, offset -> tuple(meta, bam, bai, bed12, offset) }
            // Keep model training per sample, but expose its four natural
            // stages so each expensive step can be resumed and profiled.
            MAKE_ORFRATER_TFAMS(training_inputs.map { meta, _bam, _bai, bed12, _offset -> tuple(meta, bed12) })
            FIND_ORFRATER_ORFS(MAKE_ORFRATER_TFAMS.out.tfams, fasta)
            regression_inputs = training_inputs
                .map { meta, bam, _bai, bed12, offset -> tuple(meta.id, meta, bam, bed12, offset) }
                .join(FIND_ORFRATER_ORFS.out.orfs.map { meta, orfstore, _bed12 -> tuple(meta.id, orfstore) }, by: 0)
                .map { _id, meta, bam, bed12, offset, orfstore -> tuple(meta, bam, bed12, orfstore, offset) }
            REGRESS_ORFRATER(regression_inputs)
            rate_inputs = REGRESS_ORFRATER.out.regression
                .map { meta, regression, metagene, orfstore, bed12 -> tuple(meta.id, meta, regression, metagene, orfstore, bed12) }
                .join(training_inputs.map { meta, _bam, _bai, _bed12, offset -> tuple(meta.id, offset) }, by: 0)
                .map { _id, meta, regression, metagene, orfstore, bed12, offset -> tuple(meta, regression, metagene, orfstore, bed12, offset) }
            RATE_ORFRATER(rate_inputs)
            orfrater_inputs_with_model = orfrater_inputs
                .map { meta, bam, bai, bed12 -> tuple(meta.id, meta, bam, bai, bed12) }
                .join(RATE_ORFRATER.out.model.map { meta, model -> tuple(meta.id, model) }, by: 0)
                .map { _id, meta, bam, bai, bed12, model -> tuple(meta, bam, bai, bed12, model) }
            RUN_ORFRATER(orfrater_inputs_with_model)
        }
        STANDARDISE_ORFRATER(RUN_ORFRATER.out.raw, 'orfrater', orf_gtf)
    }
    run_price = tool_selected(selected_tools, 'price') && params.partition_fai
    if (tool_selected(selected_tools, 'price')) {
        if (!params.partition_fai) {
            log.warn 'PRICE chromosome sharding requires --partition_fai; skipping PRICE while continuing other callers'
        } else {
            price_fai = file(params.partition_fai, checkIfExists: true)
            price_manifest_inputs = gn
                .map { meta, bam, bai -> tuple(meta, price_fai, bam) }
            MAKE_PRICE_CONTIG_MANIFEST(price_manifest_inputs)
            price_contigs = MAKE_PRICE_CONTIG_MANIFEST.out.manifest
                .map { meta, manifest -> tuple(meta.id, meta, manifest) }
                .flatMap { id, meta, manifest -> manifest.readLines().findAll { it && !it.startsWith('#') }.collect { line ->
                    def fields = line.split('\\t')
                    tuple(id, meta, fields[0])
                } }
            price_inputs = gn
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .combine(orf_gtf)
                .combine(fasta)
                .combine(channel.value(price_fai))
                .join(price_contigs, by: 0)
                .map { id, meta, bam, bai, gtf_file, fasta_file, fai_file, _id2, _manifest_meta, contig ->
                    tuple(meta + [shard_id: contig], bam, bai, gtf_file, fasta_file, fai_file, contig)
                }
            PREPARE_PRICE_CONTIG(price_inputs)
            price_index_inputs = PREPARE_PRICE_CONTIG.out.contig
                .map { meta, bam, bai, gtf_file, fasta_file -> tuple([id: meta.shard_id], fasta_file, gtf_file) }
            GEDI_INDEXGENOME(price_index_inputs)
            price_run_inputs = PREPARE_PRICE_CONTIG.out.contig
                .map { meta, bam, bai, gtf_file, fasta_file -> tuple(meta.shard_id, meta, bam, bai, gtf_file, fasta_file) }
                .join(GEDI_INDEXGENOME.out.index.map { meta, index -> tuple(meta.id, index) }, by: 0)
                .map { _id, meta, bam, bai, gtf_file, fasta_file, index -> tuple(meta, meta.shard_id, bam, bai, gtf_file, fasta_file, index) }
            GEDI_PRICE_SHARDED(price_run_inputs)
            STANDARDISE_PRICE(GEDI_PRICE_SHARDED.out.orfs_tsv, 'price', orf_gtf)
        }
    }
    if (run_riborf) {
        RUN_RIBORF(riborf_inputs)
        STANDARDISE_RIBORF(RUN_RIBORF.out.raw, 'riborf', orf_gtf)
    }
    if (tool_selected(selected_tools, 'ribotish')) {
        RUN_RIBOTISH(published_inputs, orf_gtf, fasta)
        STANDARDISE_RIBOTISH(RUN_RIBOTISH.out.raw, 'ribotish', orf_gtf)
    }
    if (run_ribotie) {
        if (!ribotie_gpu_enabled) {
            log.warn 'RiboTIE was selected without --ribotie_gpu true; skipping RiboTIE while continuing other callers'
            run_ribotie = false
        } else {
            ribotie_inputs = tx.ifEmpty(gn)
                .combine(orf_gtf)
                .combine(fasta)
                .map { meta, bam, bai, caller_gtf, caller_fasta -> tuple(meta, bam, bai, caller_gtf, caller_fasta) }
            PREPARE_RIBOTIE_DATA(ribotie_inputs)
            RUN_RIBOTIE(PREPARE_RIBOTIE_DATA.out.prepared)
            STANDARDISE_RIBOTIE(RUN_RIBOTIE.out.raw, 'ribotie', orf_gtf)
        }
    }

    standardized = channel.empty()
    beds = channel.empty()
    if (ribocode_enabled) {
        standardized = standardized.mix(STANDARDISE_RIBOCODE.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOCODE.out.bed12)
    }
    if (tool_selected(selected_tools, 'ribotricer')) {
        standardized = standardized.mix(STANDARDISE_RIBOTRICER.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTRICER.out.bed12)
    }
    if (tool_selected(selected_tools, 'orfquant')) {
        standardized = standardized.mix(STANDARDISE_ORFQUANT.out.standardized)
        beds = beds.mix(STANDARDISE_ORFQUANT.out.bed12)
    }
    if (run_rpbp) {
        standardized = standardized.mix(STANDARDISE_RPBP.out.standardized)
        beds = beds.mix(STANDARDISE_RPBP.out.bed12)
    }
    if (run_iribo) {
        standardized = standardized.mix(STANDARDISE_IRIBO.out.standardized)
        beds = beds.mix(STANDARDISE_IRIBO.out.bed12)
    }
    if (run_orfrater) {
        standardized = standardized.mix(STANDARDISE_ORFRATER.out.standardized)
        beds = beds.mix(STANDARDISE_ORFRATER.out.bed12)
    }
    if (run_price) {
        standardized = standardized.mix(STANDARDISE_PRICE.out.standardized)
        beds = beds.mix(STANDARDISE_PRICE.out.bed12)
    }
    if (run_riborf) {
        standardized = standardized.mix(STANDARDISE_RIBORF.out.standardized)
        beds = beds.mix(STANDARDISE_RIBORF.out.bed12)
    }
    if (tool_selected(selected_tools, 'ribotish')) {
        standardized = standardized.mix(STANDARDISE_RIBOTISH.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTISH.out.bed12)
    }
    if (run_ribotie) {
        standardized = standardized.mix(STANDARDISE_RIBOTIE.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTIE.out.bed12)
    }

    sample_tsvs = standardized
        .map { meta, _tool, path -> tuple(meta.id, path) }
        .groupTuple()
    sample_beds = beds
        .map { meta, _tool, path -> tuple(meta.id, path) }
        .groupTuple()
    sample_calls = sample_tsvs
        .join(sample_beds, by: 0)
        .map { id, tsvs, bed_files -> tuple([id: id], tsvs, bed_files) }

    if (params.run_consensus) {
        COLLECT_ORF_CALLS(sample_calls, fasta, params.min_caller_agreement)
        TRANSLON_CONSENSUS_BRIDGE(beds, gtf)
    }

    if (params.run_characterisation) {
        if (!params.run_consensus) {
            error 'run_characterisation requires run_consensus to be true'
        }
        if (!params.proteome_fasta) {
            error '--proteome_fasta is required when --run_characterisation is true'
        }
        TRANSLON_CHARACTERISATION(
            COLLECT_ORF_CALLS.out.intervals,
            gtf,
            COLLECT_ORF_CALLS.out.verdicts,
            channel.value(file(params.proteome_fasta, checkIfExists: true)),
            fasta
        )
    }
}
