#!/usr/bin/env nextflow
/*
========================================================================================
    RNASEQ PIPELINE
========================================================================================
    Short-read RNA-seq transcript assembly via STAR + StringTie2.
    Replaces the Perl/eHive StarScallopRnaseq subpipeline (using StringTie2
    instead of Scallop for the assembly step).

    Steps:
      0. (Optional) Fetch reads from ENA if no local sample sheet is provided
      1. Build STAR genome index (if not provided via --star_index)
      2. Align each RNA-seq sample with STAR (paired-end or single-end)
      3. Assemble transcripts with StringTie2 per sample
      4. Filter by coverage/length/exon count → per-sample GFF3
      5. Merge all sample GFF3 into a single non-redundant set
      6. Write output_manifest.json

    Input — provide one of:
      --sample_sheet     CSV: id,fastq_1[,fastq_2][,strandedness]
      --rnaseq_bioproject  ENA BioProject accession (reads fetched automatically)
      --rnaseq_run_accessions  Comma-separated SRR/ERR accessions

    Strandedness in sample sheet: forward | reverse | unstranded (default: unstranded)

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.fa \
        --rnaseq_bioproject PRJEB12345 \
        --outdir results/

    Run on Slurm (with local sample sheet + pre-built index):
      nextflow run . -profile slurm \
        --genome_fasta genome.fa \
        --star_index   /path/to/star_index/ \
        --sample_sheet samples.csv \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { FETCH_READS_FROM_ENA } from './modules/fetch_reads_from_ena.nf'
include { STAR_INDEX           } from './modules/star_index.nf'
include { ALIGN_AND_ASSEMBLE   } from './subworkflows/align_and_assemble.nf'
include { MERGE_RNASEQ         } from './modules/merge_rnaseq.nf'
include { WRITE_MANIFEST       } from './modules/write_manifest.nf'

params.sample_sheet              = null
params.rnaseq_bioproject         = null
params.rnaseq_run_accessions     = null
params.star_index                = null
params.annotation_gtf            = null
params.rnaseq_strandedness       = 'auto'   // auto = unstranded unless sheet says otherwise
params.max_rnaseq_runs           = 50       // cap ENA downloads
params.sjdb_overhang             = 149      // STAR sjdbOverhang (read length - 1)

def validate_params() {
    def errors = []
    if (!params.genome_fasta && !params.star_index)
        errors << "  --genome_fasta or --star_index is required"
    if (!params.sample_sheet && !params.rnaseq_bioproject && !params.rnaseq_run_accessions)
        errors << "  One of --sample_sheet, --rnaseq_bioproject, or --rnaseq_run_accessions is required"
    if (!params.outdir)
        errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

// Parse sample sheet CSV into channel of [ meta, [ reads ] ]
def parse_sample_sheet(csv_path) {
    Channel
        .fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:           row.id,
                strandedness: row.strandedness ?: 'unstranded',
            ]
            def reads = (row.fastq_2 && row.fastq_2 != '')
                ? [ file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true) ]
            [ meta, reads ]
        }
}

workflow {

    validate_params()

    // ── Step 0: Resolve sample sheet (fetch from ENA if needed) ──────────

    if (params.sample_sheet) {
        ch_sample_sheet = Channel.of(file(params.sample_sheet, checkIfExists: true))
    } else {
        // Fetch from ENA by BioProject or run accessions
        def accession = params.rnaseq_bioproject ?: params.rnaseq_run_accessions
        FETCH_READS_FROM_ENA(accession)
        ch_sample_sheet = FETCH_READS_FROM_ENA.out.sample_sheet
    }

    ch_samples = ch_sample_sheet.flatMap { sheet -> parse_sample_sheet(sheet.toString()) }

    // ── Step 1: Build or use existing STAR index ──────────────────────────

    if (params.star_index) {
        ch_star_index = Channel.of(file(params.star_index, checkIfExists: true, type: 'dir'))
    } else {
        ch_gtf = params.annotation_gtf
            ? file(params.annotation_gtf, checkIfExists: true)
            : file('NO_FILE', checkIfExists: false)
        STAR_INDEX(
            file(params.genome_fasta, checkIfExists: true),
            ch_gtf
        )
        ch_star_index = STAR_INDEX.out.index
    }

    // ── Step 2–4: Align + assemble per sample ────────────────────────────

    ALIGN_AND_ASSEMBLE(ch_samples, ch_star_index)

    // ── Step 5: Merge all sample GFF3 files ──────────────────────────────

    ch_all_gff3 = ALIGN_AND_ASSEMBLE.out.gff3
        .map { meta, gff3 -> gff3 }
        .collect()

    MERGE_RNASEQ(ch_all_gff3)

    // ── Step 6: Write manifest ────────────────────────────────────────────

    WRITE_MANIFEST(params.outdir, MERGE_RNASEQ.out.gff3)

    // ── Software versions ─────────────────────────────────────────────────

    ch_versions = Channel.empty()
        .mix(ALIGN_AND_ASSEMBLE.out.versions)
        .mix(MERGE_RNASEQ.out.versions)

    if (!params.sample_sheet) {
        ch_versions = ch_versions.mix(FETCH_READS_FROM_ENA.out.versions)
    }

    ch_versions.collectFile(
        name: 'software_versions.tsv',
        newLine: true,
        storeDir: "${params.outdir}/pipeline_info"
    )
}
