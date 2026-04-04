#!/usr/bin/env nextflow
/*
========================================================================================
    AB INITIO PIPELINE
========================================================================================
    Ab initio gene prediction with Augustus.

    Steps:
      1. Split softmasked genome FASTA into chunks
      2. Run Augustus on each chunk in parallel
      3. Parse/normalise each Augustus GFF to Ensembl-style GFF3
      4. Merge all chunks into a single GFF3 with sequential IDs
      5. Write output_manifest.json for HiveRunNextflow dataflow

    Inputs:
      --genome_fasta   Softmasked genome FASTA
      --species        Augustus species name (built-in or custom model name)
      --outdir         Output directory

    Optional:
      --augustus_config_path  Path to Augustus config/ directory (for custom models)
      --extrinsic_cfg         Extrinsic config for hints-guided prediction

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.softmasked.fa \
        --species      human \
        --outdir       results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { AUGUSTUS        } from './modules/augustus.nf'
include { PARSE_AUGUSTUS  } from './modules/parse_augustus.nf'
include { MERGE_AB_INITIO } from './modules/merge_ab_initio.nf'
include { WRITE_MANIFEST  } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.genome_fasta)
        errors << "  --genome_fasta is required (softmasked genome FASTA)"
    if (!params.species)
        errors << "  --species is required (Augustus species model name)"
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
    // SPLIT genome into chunks for parallel Augustus
    //
    ch_chunks = Channel
        .fromPath(params.genome_fasta, checkIfExists: true)
        .splitFasta(by: params.chunk_size, file: true)
        .map { fasta ->
            def idx = fasta.baseName.replaceAll(/.*_(\d+)$/, '$1')
            [ [ id: "chunk_${idx}" ], fasta ]
        }

    //
    // AUGUSTUS per chunk
    //
    AUGUSTUS(ch_chunks, params.species)

    //
    // PARSE each Augustus GFF to Ensembl-style GFF3
    //
    PARSE_AUGUSTUS(AUGUSTUS.out.gff)

    //
    // MERGE all chunks
    //
    ch_all_gff3 = PARSE_AUGUSTUS.out.gff3
        .map { meta, gff3 -> gff3 }
        .collect()

    MERGE_AB_INITIO(ch_all_gff3)

    //
    // MANIFEST for HiveRunNextflow
    //
    WRITE_MANIFEST(params.outdir, MERGE_AB_INITIO.out.gff3)

    Channel.empty()
        .mix(AUGUSTUS.out.versions)
        .mix(PARSE_AUGUSTUS.out.versions)
        .mix(MERGE_AB_INITIO.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
