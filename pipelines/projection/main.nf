#!/usr/bin/env nextflow
/*
========================================================================================
    PROJECTION PIPELINE
========================================================================================
    Project gene models from a source (reference) genome onto a target genome
    using LASTZ whole-genome alignment chains.

    Steps:
      1. Run LASTZ per target chromosome → .chain files
      2. Merge all chain files into a single sorted chain
      3. Project source GFF3 transcripts through the merged chain
      4. Write output_manifest.json for HiveRunNextflow dataflow

    Inputs:
      --source_fasta   Softmasked source (reference) genome FASTA
      --query_fasta    Target (new) genome FASTA
      --source_gff3    Source gene annotation (GFF3)
      --outdir         Output directory

    Optional (skip LASTZ if chain already exists):
      --chain          Pre-built chain file (skips LASTZ step)

    Run locally (stub):
      nextflow run . -profile local -stub \
        --source_fasta source.fa \
        --query_fasta  target.fa \
        --source_gff3  source_genes.gff3 \
        --outdir       results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { LASTZ              } from './modules/lastz.nf'
include { PROJECT_TRANSCRIPTS } from './modules/project_transcripts.nf'
include { WRITE_MANIFEST     } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.source_gff3)
        errors << "  --source_gff3 is required (source gene annotation GFF3)"
    if (!params.query_fasta && !params.chain)
        errors << "  --query_fasta or --chain is required"
    if (!params.source_fasta && !params.chain)
        errors << "  --source_fasta or --chain is required"
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
    // BUILD or USE existing chain
    //
    if (params.chain) {
        ch_chain = Channel.fromPath(params.chain, checkIfExists: true)
    } else {
        // Split source FASTA by sequence for parallel LASTZ jobs
        ch_source_seqs = Channel
            .fromPath(params.source_fasta, checkIfExists: true)
            .splitFasta(by: 1, file: true)
            .map { fasta ->
                // Use filename stem as the chromosome id
                def id = fasta.baseName.replaceAll(/\.fa(sta)?$/, '')
                [ [ id: id ], fasta ]
            }

        LASTZ(
            ch_source_seqs,
            file(params.query_fasta, checkIfExists: true)
        )

        // Merge all per-chrom chains into one
        ch_chain = LASTZ.out.chain
            .map { meta, chain -> chain }
            .collectFile(name: 'merged.chain', newLine: false)
    }

    //
    // PROJECT source GFF3 through chain
    //
    PROJECT_TRANSCRIPTS(
        ch_chain,
        file(params.source_gff3, checkIfExists: true)
    )

    //
    // MANIFEST for HiveRunNextflow
    //
    WRITE_MANIFEST(params.outdir, PROJECT_TRANSCRIPTS.out.gff3)

    Channel.empty()
        .mix(PROJECT_TRANSCRIPTS.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
