#!/usr/bin/env nextflow
/*
========================================================================================
    BEST_TARGETED PIPELINE
========================================================================================
    Targeted transcript assembly via exonerate cDNA2genome and/or protein2genome
    alignment.  Replaces the Perl/eHive BestTargetted subpipeline.

    Steps:
      1. Split cDNA FASTA into batches for parallel exonerate
      2. (Optional) Split protein FASTA into batches for parallel exonerate
      3. Align each batch with exonerate (cdna2genome or protein2genome)
      4. Filter hits by coverage and percent identity
      5. Merge all batch GFF3 results, cluster overlapping models,
         select best non-redundant set (cdna preferred over protein)
      6. Write output_manifest.json for HiveRunNextflow dataflow

    Input genome: SOFTMASKED (exonerate --softmasktarget TRUE)
    Biotype produced: best_targeted

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.fa --cdna_fasta cdnas.fa --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --genome_fasta genome.fa \
        --cdna_fasta   cdnas.fa \
        [--protein_fasta proteins.fa] \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { ALIGN_AND_FILTER    } from './subworkflows/align_and_filter.nf'
include { SELECT_BEST_TARGETED } from './modules/select_best_targeted.nf'
include { WRITE_MANIFEST       } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.genome_fasta) errors << "  --genome_fasta is required (softmasked genome)"
    if (!params.cdna_fasta && !params.protein_fasta)
        errors << "  At least one of --cdna_fasta or --protein_fasta is required"
    if (!params.outdir)       errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    ch_genome = file(params.genome_fasta, checkIfExists: true)

    ch_versions = Channel.empty()

    //
    // cDNA alignment (optional)
    //
    ch_cdna_gff3 = Channel.empty()
    if (params.cdna_fasta) {
        ch_cdna_batches = Channel
            .fromPath(params.cdna_fasta, checkIfExists: true)
            .splitFasta(by: params.cdna_batch_size, file: true)
            .map { fa ->
                def batch_id = fa.name.replaceAll(/\.fasta$|\.fa$/, '')
                [[id: batch_id], fa]
            }

        ALIGN_AND_FILTER(
            ch_cdna_batches,
            ch_genome,
            'cdna',
            params.cdna_fasta ? file(params.cdna_fasta) : null
        )
        ch_cdna_gff3 = ALIGN_AND_FILTER.out.gff3.map { meta, gff3 -> gff3 }
        ch_versions  = ch_versions.mix(ALIGN_AND_FILTER.out.versions)
    }

    //
    // Protein alignment (optional)
    //
    ch_protein_gff3 = Channel.empty()
    if (params.protein_fasta) {
        ch_protein_batches = Channel
            .fromPath(params.protein_fasta, checkIfExists: true)
            .splitFasta(by: params.protein_batch_size, file: true)
            .map { fa ->
                def batch_id = fa.name.replaceAll(/\.fasta$|\.fa$/, '')
                [[id: batch_id], fa]
            }

        ALIGN_AND_FILTER(
            ch_protein_batches,
            ch_genome,
            'protein',
            params.protein_fasta ? file(params.protein_fasta) : null
        )
        ch_protein_gff3 = ALIGN_AND_FILTER.out.gff3.map { meta, gff3 -> gff3 }
        ch_versions     = ch_versions.mix(ALIGN_AND_FILTER.out.versions)
    }

    //
    // Collect all filtered GFF3 and select best non-redundant set
    //
    ch_all_cdna_gff3    = ch_cdna_gff3.collect().ifEmpty([])
    ch_all_protein_gff3 = ch_protein_gff3.collect().ifEmpty([])

    SELECT_BEST_TARGETED(
        ch_all_cdna_gff3,
        ch_all_protein_gff3
    )
    ch_versions = ch_versions.mix(SELECT_BEST_TARGETED.out.versions)

    //
    // Manifest for HiveRunNextflow dataflow
    //
    WRITE_MANIFEST(params.outdir, SELECT_BEST_TARGETED.out.gff3)

    ch_versions.collectFile(
        name: 'software_versions.tsv',
        newLine: true,
        storeDir: "${params.outdir}/pipeline_info"
    )
}
