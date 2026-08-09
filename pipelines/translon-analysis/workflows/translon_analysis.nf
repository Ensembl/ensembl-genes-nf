include { LOAD_RIBOSEQ_OUTPUTS } from '../subworkflows/input_contract.nf'
include { RUN_RIBOCODE; RUN_RIBOTRICER; RUN_RIBOTAPER; RUN_ORFQUANT; RUN_RPBP; STANDARDISE_CALLER as STANDARDISE_RIBOCODE; STANDARDISE_CALLER as STANDARDISE_RIBOTRICER; STANDARDISE_CALLER as STANDARDISE_RIBOTAPER; STANDARDISE_CALLER as STANDARDISE_ORFQUANT; STANDARDISE_CALLER as STANDARDISE_RPBP } from '../modules/caller_run_standardise.nf'
include { COLLECT_ORF_CALLS } from '../modules/collect_orf_calls.nf'
include { TRANSLON_CONSENSUS_BRIDGE } from './translon_consensus_bridge.nf'
include { TRANSLON_CHARACTERISATION } from '../../translon-characterisation/workflows/translon_characterisation.nf'

workflow TRANSLON_ANALYSIS {
    LOAD_RIBOSEQ_OUTPUTS(
        params.riboseq_outdir,
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

    if (params.tool in ['ribocode', 'all', 'all-wave1']) {
        RUN_RIBOCODE(tx, orf_gtf, fasta)
        STANDARDISE_RIBOCODE(RUN_RIBOCODE.out.raw, 'ribocode', orf_gtf)
    }
    if (params.tool in ['ribotricer', 'all', 'all-wave1']) {
        RUN_RIBOTRICER(tx, orf_gtf, fasta)
        STANDARDISE_RIBOTRICER(RUN_RIBOTRICER.out.raw, 'ribotricer', orf_gtf)
    }
    if (params.tool in ['ribotaper', 'all', 'all-wave1']) {
        RUN_RIBOTAPER(gn, orf_gtf, fasta)
        STANDARDISE_RIBOTAPER(RUN_RIBOTAPER.out.raw, 'ribotaper', orf_gtf)
    }
    if (params.tool in ['orfquant', 'all', 'all-wave1']) {
        RUN_ORFQUANT(tx, orf_gtf, fasta)
        STANDARDISE_ORFQUANT(RUN_ORFQUANT.out.raw, 'orfquant', orf_gtf)
    }
    if (params.tool in ['rpbp', 'all', 'all-wave1']) {
        RUN_RPBP(gn, orf_gtf, fasta)
        STANDARDISE_RPBP(RUN_RPBP.out.raw, 'rpbp', orf_gtf)
    }

    standardized = channel.empty()
    beds = channel.empty()
    if (params.tool in ['ribocode', 'all', 'all-wave1']) {
        standardized = standardized.mix(STANDARDISE_RIBOCODE.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOCODE.out.bed12)
    }
    if (params.tool in ['ribotricer', 'all', 'all-wave1']) {
        standardized = standardized.mix(STANDARDISE_RIBOTRICER.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTRICER.out.bed12)
    }
    if (params.tool in ['ribotaper', 'all', 'all-wave1']) {
        standardized = standardized.mix(STANDARDISE_RIBOTAPER.out.standardized)
        beds = beds.mix(STANDARDISE_RIBOTAPER.out.bed12)
    }
    if (params.tool in ['orfquant', 'all', 'all-wave1']) {
        standardized = standardized.mix(STANDARDISE_ORFQUANT.out.standardized)
        beds = beds.mix(STANDARDISE_ORFQUANT.out.bed12)
    }
    if (params.tool in ['rpbp', 'all', 'all-wave1']) {
        standardized = standardized.mix(STANDARDISE_RPBP.out.standardized)
        beds = beds.mix(STANDARDISE_RPBP.out.bed12)
    }

    COLLECT_ORF_CALLS(
        standardized.map { _meta, _tool, path -> path }.collect(),
        beds.map { _meta, _tool, path -> path }.collect(),
        fasta,
        params.min_caller_agreement
    )

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
