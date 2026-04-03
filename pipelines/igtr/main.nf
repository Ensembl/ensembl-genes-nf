#!/usr/bin/env nextflow
/*
========================================================================================
    IGTR PIPELINE
========================================================================================
    Immunoglobulin and T-cell Receptor gene annotation via protein homology.
    Replaces the Perl/eHive IGTR_subpipeline.pm.

    Steps:
      1. Split IGTR protein FASTA (IMGT format) into parallel batches
      2. Run GenBlast per batch against the softmasked genome
      3. Filter by PID/coverage, assign ig_gene/tr_gene biotypes from IMGT headers
      4. Cluster overlapping models, select best by combined PID+Coverage score
      5. Write output_manifest.json for HiveRunNextflow dataflow

    Flat-file I/O: no Ensembl core DB dependency.
    IMGT FASTA header format: >biotype|name|source (e.g. >IG_V_gene|IGHV1-2*02|Homo sapiens)

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.softmasked.fa \
        --igtr_proteins imgt_proteins.fa \
        --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --genome_fasta genome.softmasked.fa \
        --igtr_proteins imgt_proteins.fa \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { WRITE_MANIFEST } from './modules/write_manifest.nf'
include { RUN_IGTR       } from './subworkflows/run_igtr.nf'

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

def validate_params() {
    def errors = []
    if (!params.genome_fasta)   errors << "  --genome_fasta is required (softmasked FASTA)"
    if (!params.igtr_proteins)  errors << "  --igtr_proteins is required (IMGT protein FASTA)"
    if (!params.outdir)         errors << "  --outdir is required"
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
    // GENOME — softmasked FASTA (output of repeat_masking pipeline)
    //
    ch_genome = file(params.genome_fasta, checkIfExists: true)

    //
    // IGTR PROTEINS — split into batches for parallel GenBlast
    //
    ch_protein_batches = Channel
        .fromPath(params.igtr_proteins, checkIfExists: true)
        .splitFasta(by: params.protein_batch_size, file: true)
        .map { fa ->
            def batch_id = fa.name.replaceAll(/\.fasta$|\.fa$/, '')
            [[id: batch_id], fa]
        }

    //
    // SUBWORKFLOW: GenBlast alignment + conversion + clustering
    //
    RUN_IGTR(ch_protein_batches, ch_genome)

    //
    // WRITE OUTPUT MANIFEST — required by HiveRunNextflow bridge
    //
    ch_gff3 = RUN_IGTR.out.gff3.map { meta, gff3 -> gff3 }.first()

    WRITE_MANIFEST(params.outdir, ch_gff3)

    //
    // SOFTWARE VERSIONS
    //
    RUN_IGTR.out.versions.collectFile(
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
