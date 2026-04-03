#!/usr/bin/env nextflow
/*
========================================================================================
    REFSEQ_IMPORT PIPELINE
========================================================================================
    Download and parse NCBI RefSeq GFF3 annotations into Ensembl-style GFF3.
    Replaces the Perl/eHive Refseq_import_subpipeline.

    Steps:
      1. Download RefSeq GFF3.gz from NCBI FTP using assembly accession + name
      2. Parse NCBI GFF3 → Ensembl-style GFF3 (normalise biotypes, optional
         seq-region synonym remapping of RefSeq accessions to chr names)
      3. Write output_manifest.json for HiveRunNextflow dataflow

    Run locally (stub):
      nextflow run . -profile local -stub \
        --assembly_refseq_accession GCF_000001405.40 \
        --assembly_name GRCh38.p14 \
        --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --assembly_refseq_accession GCF_000001405.40 \
        --assembly_name GRCh38.p14 \
        [--synonyms_tsv seq_region_synonyms.tsv] \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { DOWNLOAD_REFSEQ } from './modules/download_refseq.nf'
include { PARSE_REFSEQ    } from './modules/parse_refseq.nf'
include { WRITE_MANIFEST  } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.assembly_refseq_accession)
        errors << "  --assembly_refseq_accession is required (e.g. GCF_000001405.40)"
    if (!params.assembly_name)
        errors << "  --assembly_name is required (e.g. GRCh38.p14)"
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
    // DOWNLOAD RefSeq GFF3 from NCBI FTP
    //
    DOWNLOAD_REFSEQ(
        params.assembly_refseq_accession,
        params.assembly_name
    )

    //
    // PARSE NCBI GFF3 → Ensembl-style GFF3
    //
    ch_synonyms = params.synonyms_tsv
        ? file(params.synonyms_tsv, checkIfExists: true)
        : file('NO_FILE', checkIfExists: false)   // Nextflow pattern for optional file

    PARSE_REFSEQ(
        DOWNLOAD_REFSEQ.out.gff_gz,
        ch_synonyms
    )

    //
    // MANIFEST for HiveRunNextflow dataflow
    //
    WRITE_MANIFEST(params.outdir, PARSE_REFSEQ.out.gff3)

    Channel.empty()
        .mix(DOWNLOAD_REFSEQ.out.versions)
        .mix(PARSE_REFSEQ.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
