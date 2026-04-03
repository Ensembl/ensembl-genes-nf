#!/usr/bin/env nextflow
/*
========================================================================================
    LOAD_ASSEMBLY PIPELINE
========================================================================================
    Download and prepare a genome assembly from NCBI.  Replaces the relevant
    steps of the Perl/eHive LoadAssembly subpipeline (genome acquisition and
    seq-region synonym generation; DB loading is handled separately).

    Steps:
      1. Download genomic FASTA (.fna.gz) + assembly_report.txt from NCBI FTP
      2. Decompress and index genome FASTA with samtools faidx
      3. Parse assembly_report.txt → seq-region synonyms TSV + metadata JSON
      4. Write output_manifest.json for HiveRunNextflow dataflow

    Output files consumed by downstream pipelines:
      genome/  *.fna              — unmasked genome FASTA (for repeat_masking etc.)
      genome/  *.fna.fai          — samtools fai index
      genome/  *.synonyms.tsv     — RefSeq→chr_name mapping (for refseq_import)
      genome/  *.assembly_meta.json

    Run locally (stub):
      nextflow run . -profile local -stub \
        --assembly_accession GCA_000001405.29 \
        --assembly_name GRCh38.p14 \
        --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --assembly_accession GCA_000001405.29 \
        --assembly_name GRCh38.p14 \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { DOWNLOAD_ASSEMBLY    } from './modules/download_assembly.nf'
include { PREPARE_GENOME       } from './modules/prepare_genome.nf'
include { PARSE_ASSEMBLY_REPORT } from './modules/parse_assembly_report.nf'
include { WRITE_MANIFEST       } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.assembly_accession)
        errors << "  --assembly_accession is required (GCA/GCF accession, e.g. GCA_000001405.29)"
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
    // DOWNLOAD genome FASTA + assembly report from NCBI FTP
    //
    DOWNLOAD_ASSEMBLY(
        params.assembly_accession,
        params.assembly_name
    )

    //
    // INDEX genome FASTA
    //
    PREPARE_GENOME(DOWNLOAD_ASSEMBLY.out.fasta_gz)

    //
    // PARSE assembly report → synonyms + metadata
    //
    PARSE_ASSEMBLY_REPORT(DOWNLOAD_ASSEMBLY.out.report)

    //
    // MANIFEST for HiveRunNextflow dataflow
    //
    WRITE_MANIFEST(
        params.outdir,
        PREPARE_GENOME.out.fasta,
        PARSE_ASSEMBLY_REPORT.out.synonyms,
        PARSE_ASSEMBLY_REPORT.out.meta
    )

    Channel.empty()
        .mix(DOWNLOAD_ASSEMBLY.out.versions)
        .mix(PREPARE_GENOME.out.versions)
        .mix(PARSE_ASSEMBLY_REPORT.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
