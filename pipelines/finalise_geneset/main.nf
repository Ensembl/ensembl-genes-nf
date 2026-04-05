#!/usr/bin/env nextflow
/*
========================================================================================
    FINALISE_GENESET PIPELINE
========================================================================================
    Apply post-UTR-addition quality control and annotation to produce the final gene
    set.  Replaces several Ensembl eHive modules:

      HiveCleanGeneset      — removes poor-quality transcripts
      HivePseudogenes       — flags pseudogenes based on repeat coverage / intron size
      Readthrough detection — flags transcripts spanning two independent gene loci
      HiveSelenocysteineFinder — flags genes matching a selenoprotein database
      Canonical transcript selection — marks the best transcript per gene

    Inputs:
      --input_gff3              GFF3 from utr_addition (or consolidate) pipeline
      --repeat_gff3             Repeat GFF3 from repeat_masking pipeline
      --selenoprotein_fasta     (optional) FASTA of known selenoproteins
      --outdir                  Output directory

    Run locally (stub):
      nextflow run . -profile local -stub \\
        --input_gff3   genes_with_utrs.gff3 \\
        --repeat_gff3  repeats.gff3 \\
        --outdir       results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { FILTER_GENESET      } from './modules/filter_geneset.nf'
include { DETECT_PSEUDOGENES  } from './modules/detect_pseudogenes.nf'
include { DETECT_READTHROUGH  } from './modules/detect_readthrough.nf'
include { FLAG_SELENOPROTEINS } from './modules/flag_selenoproteins.nf'
include { SELECT_CANONICAL    } from './modules/select_canonical.nf'
include { WRITE_MANIFEST      } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.input_gff3)
        errors << "  --input_gff3 is required"
    if (!params.repeat_gff3)
        errors << "  --repeat_gff3 is required"
    if (!params.outdir)
        errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    //
    // Input channels
    //
    ch_input_gff3  = file(params.input_gff3,  checkIfExists: true)
    ch_repeat_gff3 = file(params.repeat_gff3, checkIfExists: true)

    // Selenoprotein FASTA is optional — pass a sentinel NO_FILE when absent so
    // FLAG_SELENOPROTEINS always runs but skips exonerate when not needed.
    ch_seleno_fasta = params.selenoprotein_fasta
        ? file(params.selenoprotein_fasta, checkIfExists: true)
        : file('NO_FILE', checkIfExists: false)

    //
    // 1. Remove poor-quality transcripts
    //
    FILTER_GENESET(
        ch_input_gff3,
        params.min_orf_aa,
        params.min_intron_size
    )

    //
    // 2. Flag pseudogenes
    //
    DETECT_PSEUDOGENES(
        FILTER_GENESET.out.gff3,
        ch_repeat_gff3,
        params.max_repeat_cds_coverage
    )

    //
    // 3. Flag readthrough transcripts
    //
    DETECT_READTHROUGH(
        DETECT_PSEUDOGENES.out.gff3,
        params.max_readthrough_gap
    )

    //
    // 4. Flag selenoproteins (always runs; no-op when NO_FILE sentinel passed)
    //
    FLAG_SELENOPROTEINS(
        DETECT_READTHROUGH.out.gff3,
        ch_seleno_fasta
    )
    ch_after_seleno = FLAG_SELENOPROTEINS.out.gff3

    //
    // 5. Select canonical transcript per gene
    //
    SELECT_CANONICAL(ch_after_seleno)

    //
    // 6. Write output manifest
    //
    WRITE_MANIFEST(params.outdir, SELECT_CANONICAL.out.gff3)

    //
    // 7. Collect software versions
    //
    Channel.empty()
        .mix(FILTER_GENESET.out.versions)
        .mix(DETECT_PSEUDOGENES.out.versions)
        .mix(DETECT_READTHROUGH.out.versions)
        .mix(FLAG_SELENOPROTEINS.out.versions)
        .mix(SELECT_CANONICAL.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
