include { LOAD_RIBOSEQ_OUTPUTS } from '../subworkflows/input_contract.nf'
include { RUN_RIBOCODE; RUN_RIBOTRICER; RUN_RIBOTAPER; RUN_ORFQUANT; RUN_RPBP; RUN_PUBLISHED_CALLER as RUN_IRIBO; RUN_PUBLISHED_CALLER as RUN_ORFRATER; RUN_PUBLISHED_CALLER as RUN_PRICE; RUN_PUBLISHED_CALLER as RUN_RIBORF; RUN_PUBLISHED_CALLER as RUN_RIBOTISH; RUN_PUBLISHED_CALLER as RUN_RIBOTIE; STANDARDISE_CALLER as STANDARDISE_RIBOCODE; STANDARDISE_CALLER as STANDARDISE_RIBOTRICER; STANDARDISE_CALLER as STANDARDISE_RIBOTAPER; STANDARDISE_CALLER as STANDARDISE_ORFQUANT; STANDARDISE_CALLER as STANDARDISE_RPBP; STANDARDISE_CALLER as STANDARDISE_IRIBO; STANDARDISE_CALLER as STANDARDISE_ORFRATER; STANDARDISE_CALLER as STANDARDISE_PRICE; STANDARDISE_CALLER as STANDARDISE_RIBORF; STANDARDISE_CALLER as STANDARDISE_RIBOTISH; STANDARDISE_CALLER as STANDARDISE_RIBOTIE } from '../modules/caller_run_standardise.nf'
include { COLLECT_ORF_CALLS } from '../modules/collect_orf_calls.nf'
include { TRANSLON_CONSENSUS_BRIDGE } from './translon_consensus_bridge.nf'
include { TRANSLON_CHARACTERISATION } from '../../translon-characterisation/workflows/translon_characterisation.nf'

def tool_enabled(selected_tools, name, wave) {
    selected_tools.contains(name) || selected_tools.contains('all') || selected_tools.contains("all-${wave}")
}

workflow TRANSLON_ANALYSIS {
    LOAD_RIBOSEQ_OUTPUTS(
        params.riboseq_outdir,
        params.samplesheet,
        params.transcriptome_bam_glob,
        params.genome_bam_glob,
        params.offsets_glob,
        params.translonscorer_glob
    )

    gtf = channel.value(file(params.gtf, checkIfExists: true))
    orf_gtf = params.canonical_gtf
        ? channel.value(file(params.canonical_gtf, checkIfExists: true))
        : gtf
    fasta = channel.value(file(params.fasta, checkIfExists: true))
    tx = LOAD_RIBOSEQ_OUTPUTS.out.transcriptome
    gn = LOAD_RIBOSEQ_OUTPUTS.out.genome
    // Published Wave 2 callers currently use one alignment contract per sample.
    // Prefer genome-space BAMs; fall back to transcriptome BAMs when necessary.
    published_inputs = gn.ifEmpty(tx)
    selected_tools = params.tools.split(',').collect { it.trim().toLowerCase() }.findAll { it }

    if (tool_enabled(selected_tools, 'ribocode', 'wave1')) {
        RUN_RIBOCODE(tx, orf_gtf, fasta)
        STANDARDISE_RIBOCODE(RUN_RIBOCODE.out.raw, 'ribocode', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'ribotricer', 'wave1')) {
        RUN_RIBOTRICER(tx, orf_gtf, fasta)
        STANDARDISE_RIBOTRICER(RUN_RIBOTRICER.out.raw, 'ribotricer', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'ribotaper', 'wave1')) {
        RUN_RIBOTAPER(gn, orf_gtf, fasta)
        STANDARDISE_RIBOTAPER(RUN_RIBOTAPER.out.raw, 'ribotaper', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'orfquant', 'wave1')) {
        RUN_ORFQUANT(tx, orf_gtf, fasta)
        STANDARDISE_ORFQUANT(RUN_ORFQUANT.out.raw, 'orfquant', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'rpbp', 'wave1')) {
        RUN_RPBP(gn, orf_gtf, fasta)
        STANDARDISE_RPBP(RUN_RPBP.out.raw, 'rpbp', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'iribo', 'wave2')) {
        RUN_IRIBO(published_inputs, orf_gtf, fasta, 'iribo')
        STANDARDISE_IRIBO(RUN_IRIBO.out.raw, 'iribo', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'orfrater', 'wave2')) {
        RUN_ORFRATER(published_inputs, orf_gtf, fasta, 'orfrater')
        STANDARDISE_ORFRATER(RUN_ORFRATER.out.raw, 'orfrater', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'price', 'wave2')) {
        RUN_PRICE(published_inputs, orf_gtf, fasta, 'price')
        STANDARDISE_PRICE(RUN_PRICE.out.raw, 'price', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'riborf', 'wave2')) {
        RUN_RIBORF(published_inputs, orf_gtf, fasta, 'riborf')
        STANDARDISE_RIBORF(RUN_RIBORF.out.raw, 'riborf', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'ribotish', 'wave2')) {
        RUN_RIBOTISH(published_inputs, orf_gtf, fasta, 'ribotish')
        STANDARDISE_RIBOTISH(RUN_RIBOTISH.out.raw, 'ribotish', orf_gtf)
    }
    if (tool_enabled(selected_tools, 'ribotie', 'wave2')) {
        RUN_RIBOTIE(published_inputs, orf_gtf, fasta, 'ribotie')
        STANDARDISE_RIBOTIE(RUN_RIBOTIE.out.raw, 'ribotie', orf_gtf)
    }

    standardized = channel.empty()
    beds = channel.empty()
    if (tool_enabled(selected_tools, 'ribocode', 'wave1')) {
        standardized = standardized.mix(STANDARDISE_RIBOCODE.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOCODE.out.bed12)
    }
    if (tool_enabled(selected_tools, 'ribotricer', 'wave1')) {
        standardized = standardized.mix(STANDARDISE_RIBOTRICER.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTRICER.out.bed12)
    }
    if (tool_enabled(selected_tools, 'ribotaper', 'wave1')) {
        standardized = standardized.mix(STANDARDISE_RIBOTAPER.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTAPER.out.bed12)
    }
    if (tool_enabled(selected_tools, 'orfquant', 'wave1')) {
        standardized = standardized.mix(STANDARDISE_ORFQUANT.out.standardized)
        beds = beds.mix(STANDARDISE_ORFQUANT.out.bed12)
    }
    if (tool_enabled(selected_tools, 'rpbp', 'wave1')) {
        standardized = standardized.mix(STANDARDISE_RPBP.out.standardized)
        beds = beds.mix(STANDARDISE_RPBP.out.bed12)
    }
    if (tool_enabled(selected_tools, 'iribo', 'wave2')) {
        standardized = standardized.mix(STANDARDISE_IRIBO.out.standardized)
        beds = beds.mix(STANDARDISE_IRIBO.out.bed12)
    }
    if (tool_enabled(selected_tools, 'orfrater', 'wave2')) {
        standardized = standardized.mix(STANDARDISE_ORFRATER.out.standardized)
        beds = beds.mix(STANDARDISE_ORFRATER.out.bed12)
    }
    if (tool_enabled(selected_tools, 'price', 'wave2')) {
        standardized = standardized.mix(STANDARDISE_PRICE.out.standardized)
        beds = beds.mix(STANDARDISE_PRICE.out.bed12)
    }
    if (tool_enabled(selected_tools, 'riborf', 'wave2')) {
        standardized = standardized.mix(STANDARDISE_RIBORF.out.standardized)
        beds = beds.mix(STANDARDISE_RIBORF.out.bed12)
    }
    if (tool_enabled(selected_tools, 'ribotish', 'wave2')) {
        standardized = standardized.mix(STANDARDISE_RIBOTISH.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTISH.out.bed12)
    }
    if (tool_enabled(selected_tools, 'ribotie', 'wave2')) {
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

    COLLECT_ORF_CALLS(sample_calls, fasta, params.min_caller_agreement)

    TRANSLON_CONSENSUS_BRIDGE(beds, gtf)

    if (params.run_characterisation) {
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
