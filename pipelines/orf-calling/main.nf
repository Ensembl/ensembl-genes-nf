#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'
validateParameters()

/*
 Inputs:
  - params.manifest OR transcriptome/genome bam dirs
  - params.gtf, params.fasta
 Behavior:
  - Loads samples as tuples [ meta, bam, bai ]
  - Dispatches to selected tool subworkflow
*/

// Load common helper to discover BAMs
workflow LOAD_SAMPLES {
    take:
    manifest
    transcriptome_bam_dir
    genome_bam_dir

    main:
    Channel
      .of(manifest)
      .filter { it }
      .map { file(it) }
      .filter { it.exists() }
      .splitCsv(header:true)
      .map { row ->
          def meta = [ id: row.sample_id as String, bam_type: (row.bam_type ?: 'transcriptome') ]
          def bam = file(row.bam, checkIfExists:true)
          def bai1 = file(bam.toString()+'.bai')
          def bai2 = file(bam.parent + '/' + bam.baseName + '.bai')
          def bai = bai1.exists() ? bai1 : (bai2.exists() ? bai2 : file('MISSING_BAI'))
          tuple(meta, bam, bai)
      }
      .set { manifest_samples }

    // Fallback: glob directories
    tm_samples = transcriptome_bam_dir ? Channel
      .fromPath("${transcriptome_bam_dir}/*.bam")
      .map { bam ->
          def bai = file(bam.toString()+'.bai'); if (!bai.exists()) bai = file(bam.parent + '/' + bam.baseName + '.bai')
          def meta = [ id: bam.baseName.replaceAll(/\.sorted$|\.toTranscriptome$/, ''), bam_type: 'transcriptome' ]
          tuple(meta, bam, bai)
      } : Channel.empty()

    gn_samples = genome_bam_dir ? Channel
      .fromPath("${genome_bam_dir}/*.bam")
      .map { bam ->
          def bai = file(bam.toString()+'.bai'); if (!bai.exists()) bai = file(bam.parent + '/' + bam.baseName + '.bai')
          def meta = [ id: bam.baseName.replaceAll(/\.sorted$/, ''), bam_type: 'genome' ]
          tuple(meta, bam, bai)
      } : Channel.empty()

    all_samples = manifest_samples.mix(tm_samples).mix(gn_samples)

    emit:
    samples = all_samples
}

include { PREP_RIBOCODE; RUN_RIBOCODE; PARSE_RIBOCODE } from './modules/ribocode.nf'
include { PREP_RIBOTRICER; RUN_RIBOTRICER; PARSE_RIBOTRICER } from './modules/ribotricer.nf'
include { PREP_RIBOTAPER; RUN_RIBOTAPER; PARSE_RIBOTAPER } from './modules/ribotaper.nf'
include { PREP_ORFQUANT; RUN_ORFQUANT; PARSE_ORFQUANT } from './modules/orfquant.nf'
include { PREP_RPBP; RUN_RPBP; PARSE_RPBP } from './modules/rpbp.nf'
// Wave 2 placeholders (real runs to be implemented)
include {
  PREP_IRIBO; RUN_IRIBO; PARSE_IRIBO;
  PREP_ORFRATER; RUN_ORFRATER; PARSE_ORFRATER;
  PREP_PRICE; RUN_PRICE; PARSE_PRICE;
  PREP_RIBORF; RUN_RIBORF; PARSE_RIBORF;
  PREP_RIBOTISH; RUN_RIBOTISH; PARSE_RIBOTISH;
  PREP_RIBOTIE; RUN_RIBOTIE; PARSE_RIBOTIE
} from './modules/wave2_placeholders.nf'
// Optional: TranslonScorer as a caller on pipeline BigWigs in future
// include { PREP_TRANSLONSCORER; RUN_TRANSLONSCORER; PARSE_TRANSLONSCORER } from './modules/translonscorer.nf'

workflow {
  def gtf_ch = Channel.value(file(params.gtf, checkIfExists:true))
  def fasta_ch = Channel.value(file(params.fasta, checkIfExists:true))
  def passdir = params.qc_pass_lengths_dir ? file(params.qc_pass_lengths_dir) : null

  LOAD_SAMPLES(params.manifest, params.transcriptome_bam_dir, params.genome_bam_dir)

  def tool = params.tool ?: 'all'

  // Split samples by bam_type to route to correct tools
  def samples_tx = LOAD_SAMPLES.out.samples.filter { it[0].bam_type == 'transcriptome' }
  def samples_gn = LOAD_SAMPLES.out.samples.filter { it[0].bam_type == 'genome' }

  // RiboCode (transcriptome)
  if (tool in ['ribocode','all','all-wave1']) {
    PREP_RIBOCODE(samples_tx, gtf_ch, fasta_ch)
    RUN_RIBOCODE(PREP_RIBOCODE.out.prepared)
    PARSE_RIBOCODE(RUN_RIBOCODE.out.raw)
  }
  // ribotricer (transcriptome)
  if (tool in ['ribotricer','all','all-wave1']) {
    PREP_RIBOTRICER(samples_tx, gtf_ch, fasta_ch)
    RUN_RIBOTRICER(PREP_RIBOTRICER.out.prepared)
    PARSE_RIBOTRICER(RUN_RIBOTRICER.out.raw)
  }
  // RiboTaper (genome)
  if (tool in ['ribotaper','all','all-wave1']) {
    PREP_RIBOTAPER(samples_gn, gtf_ch, fasta_ch)
    RUN_RIBOTAPER(PREP_RIBOTAPER.out.prepared)
    PARSE_RIBOTAPER(RUN_RIBOTAPER.out.raw)
  }
  // ORFquant (transcriptome)
  if (tool in ['orfquant','all','all-wave1']) {
    PREP_ORFQUANT(samples_tx, gtf_ch, fasta_ch)
    RUN_ORFQUANT(PREP_ORFQUANT.out.prepared)
    PARSE_ORFQUANT(RUN_ORFQUANT.out.raw)
  }
  // Rp-Bp (genome)
  if (tool in ['rpbp','all','all-wave1']) {
    PREP_RPBP(samples_gn, gtf_ch, fasta_ch)
    RUN_RPBP(PREP_RPBP.out.prepared)
    PARSE_RPBP(RUN_RPBP.out.raw)
  }

  // Wave 2 (placeholders for now)
  if (tool in ['iribo','all','all-wave2']) {
    PREP_IRIBO(samples_gn.mix(samples_tx), gtf_ch, fasta_ch)
    RUN_IRIBO(PREP_IRIBO.out.prepared)
    PARSE_IRIBO(RUN_IRIBO.out.raw)
  }
  if (tool in ['orfrater','all','all-wave2']) {
    PREP_ORFRATER(samples_gn.mix(samples_tx), gtf_ch, fasta_ch)
    RUN_ORFRATER(PREP_ORFRATER.out.prepared)
    PARSE_ORFRATER(RUN_ORFRATER.out.raw)
  }
  if (tool in ['price','all','all-wave2']) {
    PREP_PRICE(samples_gn.mix(samples_tx), gtf_ch, fasta_ch)
    RUN_PRICE(PREP_PRICE.out.prepared)
    PARSE_PRICE(RUN_PRICE.out.raw)
  }
  if (tool in ['riborf','all','all-wave2']) {
    PREP_RIBORF(samples_gn.mix(samples_tx), gtf_ch, fasta_ch)
    RUN_RIBORF(PREP_RIBORF.out.prepared)
    PARSE_RIBORF(RUN_RIBORF.out.raw)
  }
  if (tool in ['ribotish','all','all-wave2']) {
    PREP_RIBOTISH(samples_gn.mix(samples_tx), gtf_ch, fasta_ch)
    RUN_RIBOTISH(PREP_RIBOTISH.out.prepared)
    PARSE_RIBOTISH(RUN_RIBOTISH.out.raw)
  }
  if (tool in ['ribotie','all','all-wave2']) {
    PREP_RIBOTIE(samples_gn.mix(samples_tx), gtf_ch, fasta_ch)
    RUN_RIBOTIE(PREP_RIBOTIE.out.prepared)
    PARSE_RIBOTIE(RUN_RIBOTIE.out.raw)
  }
}
