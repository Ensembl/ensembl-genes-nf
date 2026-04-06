#!/usr/bin/env nextflow
/*
========================================================================================
    GFF3_TO_CORE PIPELINE
========================================================================================
    Loads a final-geneset GFF3 file into an Ensembl core MySQL database schema,
    then assigns Ensembl stable IDs (ENSG/ENST/ENSE/ENSP) to all features.

    This is the final stage in the Nextflow genebuild pipeline, replacing the
    eHive-based core DB loading that uses the Perl Ensembl API.

    Stages:
      1. LOAD_GFF3_TO_CORE   — parse GFF3 → INSERT into gene/transcript/exon/
                               translation/exon_transcript tables
      2. ASSIGN_STABLE_IDS   — assign ENS*G/T/E/P stable IDs sequentially
      3. WRITE_MANIFEST       — record completion + DB name in output_manifest.json

    Prerequisites:
      - An empty (or pre-initialized) Ensembl core schema DB must exist.
        Create it with:
          mysql -u ensadmin -p -e "CREATE DATABASE homo_sapiens_core_109_38;"
          mysql -u ensadmin -p homo_sapiens_core_109_38 < ensembl-core-schema.sql
      - The genome .fai index (from samtools faidx) for seq_region loading
      - The synonyms TSV from load_assembly (for seq-region name mapping)

    Required parameters:
      --input_gff3   Final geneset GFF3
      --db_host      MySQL host
      --db_user      MySQL user (needs INSERT/UPDATE privileges)
      --db_name      Database name
      --assembly     Assembly version string (e.g. GRCh38)
      --outdir       Output directory for stats and manifest

    Optional parameters:
      --genome_fai          samtools .fai index path (for seq_region loading)
      --synonyms_tsv        Synonyms TSV from load_assembly
      --stable_id_prefix    '' for human, 'GAL' for chicken, etc.
      --species_name        Ensembl production name (e.g. homo_sapiens)

    Run locally (stub):
      nextflow run . -profile local -stub \\
        --input_gff3 final.gff3 \\
        --db_host localhost --db_user ensadmin --db_name hs_core \\
        --assembly GRCh38 --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { LOAD_GFF3_TO_CORE  } from './modules/load_genes.nf'
include { ASSIGN_STABLE_IDS  } from './modules/assign_stable_ids.nf'
include { WRITE_MANIFEST     } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.input_gff3) errors << "  --input_gff3 is required"
    if (!params.db_host)    errors << "  --db_host is required"
    if (!params.db_user)    errors << "  --db_user is required"
    if (!params.db_name)    errors << "  --db_name is required"
    if (!params.assembly)   errors << "  --assembly is required"
    if (!params.outdir)     errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    ch_gff3        = file(params.input_gff3, checkIfExists: true)
    ch_genome_fai  = params.genome_fai   ? file(params.genome_fai,   checkIfExists: true)
                                         : file("NO_FILE_FAI")
    ch_synonyms    = params.synonyms_tsv ? file(params.synonyms_tsv, checkIfExists: true)
                                         : file("NO_FILE_SYN")

    // Stage 1: Load genes
    LOAD_GFF3_TO_CORE(ch_gff3, ch_genome_fai, ch_synonyms)

    // Stage 2: Assign stable IDs (depends on Stage 1 completing first)
    ASSIGN_STABLE_IDS(LOAD_GFF3_TO_CORE.out.stats)

    // Stage 3: Write manifest
    WRITE_MANIFEST(
        params.outdir,
        LOAD_GFF3_TO_CORE.out.stats,
        ASSIGN_STABLE_IDS.out.stats,
    )

    Channel.empty()
        .mix(LOAD_GFF3_TO_CORE.out.versions)
        .mix(ASSIGN_STABLE_IDS.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
