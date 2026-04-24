#!/usr/bin/env nextflow
/*
========================================================================================
    CONSOLIDATE PIPELINE
========================================================================================
    Merge gene models from multiple annotation sources into a single geneset
    using a layer-annotation priority strategy.

    Priority layers (lower integer = higher quality evidence):
      0  long_read        IsoSeq / long-read direct evidence
      1  best_targeted    cDNA/protein best-targeted
      2  rnaseq           short-read assembled
      3  projection       cross-genome projection
      4  refseq           RefSeq/external annotation
      5  ab_initio        Augustus ab initio
      6  genblast         genblast homology

    A lower-priority transcript is suppressed if it overlaps a retained
    higher-priority transcript; otherwise it is kept (fills in the gaps).

    Inputs:
      --gff3_dir          Directory containing GFF3 files to consolidate
                          OR individual files via --gff3_files
      --layer_priorities  JSON map of filename-pattern → priority integer
                          e.g. '{"rnaseq":2,"ab_initio":5,"refseq":4}'
      --outdir            Output directory

    Run locally (stub):
      nextflow run . -profile local -stub \
        --gff3_dir      annotation_sources/ \
        --layer_priorities '{"rnaseq":2,"ab_initio":5}' \
        --outdir        results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { CONSOLIDATE_GENES } from './modules/consolidate_genes.nf'
include { WRITE_MANIFEST    } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.gff3_dir && !params.gff3_files)
        errors << "  --gff3_dir or --gff3_files is required"
    if (!params.layer_priorities)
        errors << "  --layer_priorities is required (JSON: {\"pattern\": priority})"
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
    // Collect GFF3 files
    // --gff3_files accepts:
    //   • a single path
    //   • a comma-separated list of paths (used when master pipeline collects outputs)
    //   • a glob pattern
    // --gff3_dir triggers a recursive glob over the directory.
    //
    if (params.gff3_files) {
        def gff3_list = params.gff3_files instanceof String
            ? params.gff3_files.split(',').collect { it.trim() }.findAll { it }
            : [params.gff3_files.toString()]

        ch_gff3 = Channel
            .fromList(gff3_list)
            .map { p -> file(p, checkIfExists: true) }
            .collect()
    } else {
        // Recursive glob — each annotation pipeline writes its GFF3 one or
        // two levels below outdir (e.g. outdir/rnaseq/rnaseq.gff3).
        // Repeat GFF3s are NOT published to the outdir so they won't appear.
        ch_gff3 = Channel
            .fromPath("${params.gff3_dir}/**/*.gff3", checkIfExists: true)
            .collect()
    }

    //
    // CONSOLIDATE
    //
    CONSOLIDATE_GENES(ch_gff3)

    //
    // MANIFEST
    //
    WRITE_MANIFEST(params.outdir, CONSOLIDATE_GENES.out.gff3)

    Channel.empty()
        .mix(CONSOLIDATE_GENES.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
