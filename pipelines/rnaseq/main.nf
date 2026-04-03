#!/usr/bin/env nextflow
/*
========================================================================================
    RNASEQ PIPELINE
========================================================================================
    Short-read RNA-seq transcript assembly via STAR + StringTie2.
    Replaces the Perl/eHive StarScallopRnaseq subpipeline (using StringTie2
    instead of Scallop for the assembly step).

    Steps:
      1. Build STAR genome index (if not provided)
      2. Align each RNA-seq sample with STAR (paired-end or single-end)
      3. Assemble transcripts with StringTie2 per sample
      4. Filter by coverage/length/exon count → per-sample GFF3
      5. Merge all sample GFF3 into a single non-redundant set
      6. Write output_manifest.json for HiveRunNextflow dataflow

    Input sample sheet (CSV): id,fastq_1[,fastq_2][,strandedness]
      strandedness: forward | reverse | unstranded (default: unstranded)

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.fa \
        --sample_sheet samples.csv \
        --outdir results/

    Run on Slurm (with pre-built index):
      nextflow run . -profile slurm \
        --genome_fasta genome.fa \
        --star_index   /path/to/star_index/ \
        --sample_sheet samples.csv \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { STAR_INDEX        } from './modules/star_index.nf'
include { ALIGN_AND_ASSEMBLE } from './subworkflows/align_and_assemble.nf'
include { MERGE_RNASEQ      } from './modules/merge_rnaseq.nf'
include { WRITE_MANIFEST    } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.genome_fasta && !params.star_index)
        errors << "  --genome_fasta or --star_index is required"
    if (!params.sample_sheet)
        errors << "  --sample_sheet is required (CSV: id,fastq_1[,fastq_2][,strandedness])"
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
                id:            row.id,
                strandedness:  row.strandedness ?: 'unstranded',
            ]
            def reads = row.fastq_2
                ? [ file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true) ]
            [ meta, reads ]
        }
}

workflow {

    validate_params()

    //
    // BUILD or USE existing STAR index
    //
    if (params.star_index) {
        ch_star_index = file(params.star_index, checkIfExists: true, type: 'dir')
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

    //
    // ALIGN + ASSEMBLE per sample
    //
    ch_samples = parse_sample_sheet(params.sample_sheet)

    ALIGN_AND_ASSEMBLE(ch_samples, ch_star_index)

    //
    // MERGE all sample GFF3 files
    //
    ch_all_gff3 = ALIGN_AND_ASSEMBLE.out.gff3
        .map { meta, gff3 -> gff3 }
        .collect()

    MERGE_RNASEQ(ch_all_gff3)

    //
    // MANIFEST for HiveRunNextflow
    //
    WRITE_MANIFEST(params.outdir, MERGE_RNASEQ.out.gff3)

    Channel.empty()
        .mix(ALIGN_AND_ASSEMBLE.out.versions)
        .mix(MERGE_RNASEQ.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
