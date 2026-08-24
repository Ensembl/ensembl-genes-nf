include { LOAD_RIBOSEQ_OUTPUTS } from '../subworkflows/input_contract.nf'
include { MERGE_RIBO_INPUTS; INFLATE_UNIQUE_BAM } from '../modules/merge_bams.nf'
include { PREPARE_CALLER_INPUTS } from '../subworkflows/caller_inputs.nf'
include { RUN_RIBOCODE; RUN_RIBOTRICER; RUN_ORFQUANT; RUN_RPBP; RUN_IRIBO; TRAIN_ORFRATER; RUN_ORFRATER; RUN_RIBORF; RUN_RIBOTISH; RUN_RIBOTIE; STANDARDISE_CALLER as STANDARDISE_RIBOCODE; STANDARDISE_CALLER as STANDARDISE_RIBOTRICER; STANDARDISE_CALLER as STANDARDISE_ORFQUANT; STANDARDISE_CALLER as STANDARDISE_RPBP; STANDARDISE_CALLER as STANDARDISE_IRIBO; STANDARDISE_CALLER as STANDARDISE_ORFRATER; STANDARDISE_CALLER as STANDARDISE_PRICE; STANDARDISE_CALLER as STANDARDISE_RIBORF; STANDARDISE_CALLER as STANDARDISE_RIBOTISH; STANDARDISE_CALLER as STANDARDISE_RIBOTIE } from '../modules/caller_run_standardise.nf'
include { GEDI_INDEXGENOME } from '../../../modules/nf-core/gedi/indexgenome/main.nf'
include { GEDI_PRICE } from '../../../modules/nf-core/gedi/price/main.nf'
include { COLLECT_ORF_CALLS } from '../modules/collect_orf_calls.nf'
include { TRANSLON_CONSENSUS_BRIDGE } from './translon_consensus_bridge.nf'
include { TRANSLON_CHARACTERISATION } from '../../translon-characterisation/workflows/translon_characterisation.nf'
include { RIBOMETRIC } from '../../riboseq/modules/ribometric.nf'

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
        .filter { meta, bam, bai -> meta.bam_type == 'transcriptome' }
    gn = INFLATE_UNIQUE_BAM.out.inflated
        .filter { meta, bam, bai -> meta.bam_type == 'genome' }
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
    selected_tools = params.tools.split(',').collect { it.trim().toLowerCase() }.findAll { it }
    // The upstream transcriptome BAM is already the native input for
    // Ribotricer, RiboTIE and most learned callers. Prepare legacy models
    // only for callers that require genePred/BED/SAM representations.
    transcript_models = channel.empty()
    riborf_inputs = channel.empty()
    if (tool_selected(selected_tools, 'orfrater') || tool_selected(selected_tools, 'riborf')) {
        if (tool_selected(selected_tools, 'riborf') && !params.samplesheet) {
            error 'RibORF requires a samplesheet with an offsets column so RiboSeq QC offsets can be converted for offsetCorrect.pl'
        }
        // Do not call ifEmpty here. In merged mode this channel is produced by
        // RiboMetric, so checking it during workflow construction races the
        // upstream process and can falsely report that offsets are missing.
        caller_offsets = tool_selected(selected_tools, 'riborf') ? offsets : channel.empty()
        PREPARE_CALLER_INPUTS(tx, orf_gtf, caller_offsets)
        transcript_models = PREPARE_CALLER_INPUTS.out.transcript_models
        riborf_inputs = PREPARE_CALLER_INPUTS.out.riborf
    }
    // Published callers currently use one alignment contract per sample.
    // Prefer genome-space BAMs; fall back to transcriptome BAMs when necessary.
    published_inputs = gn.ifEmpty(tx)

    if (tool_selected(selected_tools, 'ribocode')) {
        RUN_RIBOCODE(tx, orf_gtf, fasta)
        STANDARDISE_RIBOCODE(RUN_RIBOCODE.out.raw, 'ribocode', orf_gtf)
    }
    if (tool_selected(selected_tools, 'ribotricer')) {
        RUN_RIBOTRICER(tx, orf_gtf, fasta)
        STANDARDISE_RIBOTRICER(RUN_RIBOTRICER.out.raw, 'ribotricer', orf_gtf)
    }
    if (tool_selected(selected_tools, 'orfquant')) {
        orfquant_alignments = gn.ifEmpty(tx)
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .join(offsets.map { meta, offset -> tuple(meta.id, offset) }, by: 0)
            .map { id, meta, bam, bai, offset -> tuple(meta, bam, bai, offset) }
        RUN_ORFQUANT(orfquant_alignments, orf_gtf, fasta)
        STANDARDISE_ORFQUANT(RUN_ORFQUANT.out.raw, 'orfquant', orf_gtf)
    }
    if (tool_selected(selected_tools, 'rpbp') && !params.skip_fastq_tools) {
        if (!params.ribosomal_fasta || !params.adapter_fasta) {
            error 'Rp-Bp was selected but --ribosomal_fasta and --adapter_fasta were not provided'
        }
        RUN_RPBP(ribo_fastq.ifEmpty { error 'Rp-Bp was selected but no ribo_fastq is present in the samplesheet' }, orf_gtf, fasta, ribosomal_fasta, adapter_fasta)
        STANDARDISE_RPBP(RUN_RPBP.out.raw, 'rpbp', orf_gtf)
    }
    if (tool_selected(selected_tools, 'iribo')) {
        RUN_IRIBO(published_inputs, orf_gtf, fasta)
        STANDARDISE_IRIBO(RUN_IRIBO.out.raw, 'iribo', orf_gtf)
    }
    if (tool_selected(selected_tools, 'orfrater')) {
        orfrater_inputs = transcript_models
            .map { meta, bam, bai, _genepred, bed12 -> tuple(meta, bam, bai, bed12) }
        if (params.orfrater_model) {
            orfrater_model_path = file(params.orfrater_model, checkIfExists: true)
            required_orfrater_files = ['orfratings.h5', 'metagene.txt', 'offsets.txt']
            missing_orfrater_files = required_orfrater_files.findAll { name -> !file("${orfrater_model_path}/${name}").exists() }
            if (missing_orfrater_files) {
                error "ORF-RATER model directory is missing: ${missing_orfrater_files.join(', ')} (${orfrater_model_path})"
            }
            orfrater_inputs_with_model = orfrater_inputs
                .map { meta, bam, bai, bed12 -> tuple(meta, bam, bai, bed12, orfrater_model_path) }
            RUN_ORFRATER(orfrater_inputs_with_model, orf_gtf, fasta)
        } else {
            training_inputs = transcript_models
                .map { meta, bam, bai, _genepred, bed12 -> tuple(meta.id, meta, bam, bai, bed12) }
                .join(offsets.map { meta, offset -> tuple(meta.id, offset) }, by: 0)
                .map { id, meta, bam, bai, bed12, offset -> tuple(meta, bam, bai, bed12, offset) }
            TRAIN_ORFRATER(training_inputs, fasta)
            orfrater_inputs_with_model = orfrater_inputs
                .map { meta, bam, bai, bed12 -> tuple(meta.id, meta, bam, bai, bed12) }
                .join(TRAIN_ORFRATER.out.model.map { meta, model -> tuple(meta.id, model) }, by: 0)
                .map { id, meta, bam, bai, bed12, model -> tuple(meta, bam, bai, bed12, model) }
            RUN_ORFRATER(orfrater_inputs_with_model, orf_gtf, fasta)
        }
        STANDARDISE_ORFRATER(RUN_ORFRATER.out.raw, 'orfrater', orf_gtf)
    }
    if (tool_selected(selected_tools, 'price')) {
        price_index_input = fasta
            .combine(orf_gtf)
            .map { ref_fasta, ref_gtf -> tuple([id: params.gedi_reference_id ?: 'reference'], ref_fasta, ref_gtf) }
        GEDI_INDEXGENOME(price_index_input)
        price_bam_input = gn
            .map { meta, bam, bai -> tuple('price_cohort', bam, bai) }
            .groupTuple()
            .map { cohort, bams, bais -> tuple([id: cohort], bams, bais) }
        GEDI_PRICE(price_bam_input, GEDI_INDEXGENOME.out.index)
        STANDARDISE_PRICE(GEDI_PRICE.out.orfs_tsv, 'price', orf_gtf)
    }
    if (tool_selected(selected_tools, 'riborf')) {
        RUN_RIBORF(riborf_inputs, orf_gtf, fasta)
        STANDARDISE_RIBORF(RUN_RIBORF.out.raw, 'riborf', orf_gtf)
    }
    if (tool_selected(selected_tools, 'ribotish')) {
        RUN_RIBOTISH(published_inputs, orf_gtf, fasta)
        STANDARDISE_RIBOTISH(RUN_RIBOTISH.out.raw, 'ribotish', orf_gtf)
    }
    if (tool_selected(selected_tools, 'ribotie')) {
        if (!ribotie_gpu_enabled) {
            error 'RiboTIE requires --ribotie_gpu true and a CUDA-capable container'
        }
        RUN_RIBOTIE(published_inputs, orf_gtf, fasta)
        STANDARDISE_RIBOTIE(RUN_RIBOTIE.out.raw, 'ribotie', orf_gtf)
    }

    standardized = channel.empty()
    beds = channel.empty()
    if (tool_selected(selected_tools, 'ribocode')) {
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
    if (tool_selected(selected_tools, 'rpbp') && !params.skip_fastq_tools) {
        standardized = standardized.mix(STANDARDISE_RPBP.out.standardized)
        beds = beds.mix(STANDARDISE_RPBP.out.bed12)
    }
    if (tool_selected(selected_tools, 'iribo')) {
        standardized = standardized.mix(STANDARDISE_IRIBO.out.standardized)
        beds = beds.mix(STANDARDISE_IRIBO.out.bed12)
    }
    if (tool_selected(selected_tools, 'orfrater')) {
        standardized = standardized.mix(STANDARDISE_ORFRATER.out.standardized)
        beds = beds.mix(STANDARDISE_ORFRATER.out.bed12)
    }
    if (tool_selected(selected_tools, 'price')) {
        standardized = standardized.mix(STANDARDISE_PRICE.out.standardized)
        beds = beds.mix(STANDARDISE_PRICE.out.bed12)
    }
    if (tool_selected(selected_tools, 'riborf')) {
        standardized = standardized.mix(STANDARDISE_RIBORF.out.standardized)
        beds = beds.mix(STANDARDISE_RIBORF.out.bed12)
    }
    if (tool_selected(selected_tools, 'ribotish')) {
        standardized = standardized.mix(STANDARDISE_RIBOTISH.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTISH.out.bed12)
    }
    if (tool_selected(selected_tools, 'ribotie')) {
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
