include { LOAD_RIBOSEQ_OUTPUTS } from '../subworkflows/input_contract.nf'
include { PREP_RIBOCODE; RUN_RIBOCODE; PARSE_RIBOCODE } from '../../orf-calling/modules/ribocode.nf'
include { PREP_RIBOTRICER; RUN_RIBOTRICER; PARSE_RIBOTRICER } from '../../orf-calling/modules/ribotricer.nf'
include { PREP_RIBOTAPER; RUN_RIBOTAPER; PARSE_RIBOTAPER } from '../../orf-calling/modules/ribotaper.nf'
include { PREP_ORFQUANT; RUN_ORFQUANT; PARSE_ORFQUANT } from '../../orf-calling/modules/orfquant.nf'
include { PREP_RPBP; RUN_RPBP; PARSE_RPBP } from '../../orf-calling/modules/rpbp.nf'
include { COLLECT_ORF_CALLS } from '../modules/collect_orf_calls.nf'
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
    fasta = channel.value(file(params.fasta, checkIfExists: true))
    tx = LOAD_RIBOSEQ_OUTPUTS.out.transcriptome
    gn = LOAD_RIBOSEQ_OUTPUTS.out.genome

    if (params.tool in ['ribocode', 'all', 'all-wave1']) {
        PREP_RIBOCODE(tx, gtf, fasta)
        RUN_RIBOCODE(PREP_RIBOCODE.out.prepared)
        PARSE_RIBOCODE(RUN_RIBOCODE.out.raw)
    }
    if (params.tool in ['ribotricer', 'all', 'all-wave1']) {
        PREP_RIBOTRICER(tx, gtf, fasta)
        RUN_RIBOTRICER(PREP_RIBOTRICER.out.prepared)
        PARSE_RIBOTRICER(RUN_RIBOTRICER.out.raw)
    }
    if (params.tool in ['ribotaper', 'all', 'all-wave1']) {
        PREP_RIBOTAPER(gn, gtf, fasta)
        RUN_RIBOTAPER(PREP_RIBOTAPER.out.prepared)
        PARSE_RIBOTAPER(RUN_RIBOTAPER.out.raw)
    }
    if (params.tool in ['orfquant', 'all', 'all-wave1']) {
        PREP_ORFQUANT(tx, gtf, fasta)
        RUN_ORFQUANT(PREP_ORFQUANT.out.prepared)
        PARSE_ORFQUANT(RUN_ORFQUANT.out.raw)
    }
    if (params.tool in ['rpbp', 'all', 'all-wave1']) {
        PREP_RPBP(gn, gtf, fasta)
        RUN_RPBP(PREP_RPBP.out.prepared)
        PARSE_RPBP(RUN_RPBP.out.raw)
    }

    standardized = channel.empty()
    beds = channel.empty()
    if (params.tool in ['ribocode', 'all', 'all-wave1']) {
        standardized = standardized.mix(PARSE_RIBOCODE.out.standardized)
        beds = beds.mix(PARSE_RIBOCODE.out.bed12)
    }
    if (params.tool in ['ribotricer', 'all', 'all-wave1']) {
        standardized = standardized.mix(PARSE_RIBOTRICER.out.standardized)
        beds = beds.mix(PARSE_RIBOTRICER.out.bed12)
    }
    if (params.tool in ['ribotaper', 'all', 'all-wave1']) {
        standardized = standardized.mix(PARSE_RIBOTAPER.out.standardized)
        beds = beds.mix(PARSE_RIBOTAPER.out.bed12)
    }
    if (params.tool in ['orfquant', 'all', 'all-wave1']) {
        standardized = standardized.mix(PARSE_ORFQUANT.out.standardized)
        beds = beds.mix(PARSE_ORFQUANT.out.bed12)
    }
    if (params.tool in ['rpbp', 'all', 'all-wave1']) {
        standardized = standardized.mix(PARSE_RPBP.out.standardized)
        beds = beds.mix(PARSE_RPBP.out.bed12)
    }

    COLLECT_ORF_CALLS(
        standardized.map { _meta, path -> path }.collect(),
        beds.map { _meta, path -> path }.collect(),
        fasta,
        params.min_caller_agreement
    )

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
