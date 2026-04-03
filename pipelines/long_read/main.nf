#!/usr/bin/env nextflow
/*
========================================================================================
    LONG_READ PIPELINE
========================================================================================
    Replaces the Perl/eHive long_read subpipeline.
    Aligns long reads (ONT / PacBio) with minimap2, collapses transcript models,
    and classifies by UniProt protein support. Outputs GFF3 + output_manifest.json.

    Flat-file I/O: no Ensembl core DB dependency.
    output_manifest.json is read by HiveRunNextflow to fan out downstream eHive jobs.

    Run locally (stub, no containers):
      nextflow run . -profile local -stub \
        --sample_sheet samples.tsv --genome_fasta genome.fa \
        --protein_db /path/to/uniprot_db --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --sample_sheet samples.tsv --genome_fasta genome.fa \
        --protein_db /path/to/uniprot_db --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { MINIMAP2_INDEX } from './modules/minimap2_index.nf'
include { WRITE_MANIFEST } from './modules/write_manifest.nf'
include { ALIGN          } from './subworkflows/align.nf'
include { COLLAPSE       } from './subworkflows/collapse.nf'
include { CLASSIFY       } from './subworkflows/classify.nf'

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

def validate_params() {
    def errors = []
    if (!params.sample_sheet) errors << "  --sample_sheet is required"
    if (!params.genome_fasta)  errors << "  --genome_fasta is required"
    if (!params.protein_db)    errors << "  --protein_db is required"
    if (!params.outdir)        errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    validate_params()

    //
    // INPUT: Sample sheet
    // TSV with columns: sample_name, fastq_file, instrument_platform (ONT | PACBIO_SMRT)
    //
    ch_samples = Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id:                  row.sample_name,
                instrument_platform: row.instrument_platform ?: 'ONT',
                minimap2_preset:     row.instrument_platform == 'PACBIO_SMRT'
                                     ? 'splice:hq' : params.minimap2_preset
            ]
            def fastq = file(row.fastq_file, checkIfExists: !params.skip_download)
            [meta, fastq]
        }

    //
    // REFERENCE GENOME
    //
    ch_genome = Channel.of([[id: 'genome'], file(params.genome_fasta, checkIfExists: true)])

    //
    // MINIMAP2 INDEX — build or load pre-built index
    //
    if (params.genome_index) {
        ch_index = Channel.of([[id: 'genome'], file(params.genome_index, checkIfExists: true)])
    } else {
        MINIMAP2_INDEX(ch_genome)
        ch_index = MINIMAP2_INDEX.out.index
    }

    //
    // SUBWORKFLOW: Align reads to genome
    //
    ALIGN(ch_samples, ch_index)

    //
    // SUBWORKFLOW: Collapse overlapping alignments into transcript models
    //
    COLLAPSE(
        ALIGN.out.bam_bai,
        params.collapse_min_overlap,
        params.max_intron_size
    )

    //
    // SUBWORKFLOW: Classify by protein support (blastp against UniProt)
    //
    ch_protein_db = Channel.of([[id: 'uniprot'], file(params.protein_db, type: 'dir', checkIfExists: true)])

    CLASSIFY(
        COLLAPSE.out.gff3,
        ch_protein_db
    )

    //
    // WRITE OUTPUT MANIFEST — required by HiveRunNextflow bridge
    //
    ch_gff3_all = CLASSIFY.out.gff3.map { meta, gff3 -> gff3 }.collect()
    ch_bam_all  = ALIGN.out.bam.map     { meta, bam  -> bam  }.collect()

    WRITE_MANIFEST(
        params.outdir,
        ch_gff3_all,
        ch_bam_all
    )

    //
    // SOFTWARE VERSIONS
    //
    ch_versions = Channel.empty()
        .mix(ALIGN.out.versions)
        .mix(COLLAPSE.out.versions)
        .mix(CLASSIFY.out.versions)

    ch_versions.collectFile(
        name: 'software_versions.tsv',
        newLine: true,
        storeDir: "${params.outdir}/pipeline_info"
    )
}

/*
========================================================================================
    THE END
========================================================================================
*/
