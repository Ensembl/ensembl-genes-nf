#!/usr/bin/env nextflow
/*
========================================================================================
    GENBLAST_HOMOLOGY PIPELINE
========================================================================================
    Protein homology gene annotation via GenBlast against UniProt proteins.
    Replaces the Perl/eHive GenblastHomology subpipeline.

    Steps:
      1. Split UniProt FASTA into parallel batches
      2. Run GenBlast per batch against the softmasked genome
      3. Classify models into genblast_1..7 tiers (by PID + coverage)
      4. Merge and remove redundant models by biotype priority
      5. Write output_manifest.json for HiveRunNextflow dataflow

    Flat-file I/O: no Ensembl core DB dependency.
    Use --uniprot_set to document which clade protein set was used (informational).

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.softmasked.fa \
        --uniprot_fasta uniprot_mammals.fa \
        --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { WRITE_MANIFEST } from './modules/write_manifest.nf'
include { RUN_HOMOLOGY   } from './subworkflows/run_homology.nf'

def validate_params() {
    def errors = []
    if (!params.genome_fasta)   errors << "  --genome_fasta is required"
    if (!params.uniprot_fasta)  errors << "  --uniprot_fasta is required"
    if (!params.outdir)         errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    ch_genome = file(params.genome_fasta, checkIfExists: true)

    // Split UniProt FASTA into parallel batches
    ch_protein_batches = Channel
        .fromPath(params.uniprot_fasta, checkIfExists: true)
        .splitFasta(by: params.protein_batch_size, file: true)
        .map { fa ->
            def batch_id = fa.name.replaceAll(/\.fasta$|\.fa$/, '')
            [[id: batch_id], fa]
        }

    RUN_HOMOLOGY(ch_protein_batches, ch_genome)

    ch_gff3 = RUN_HOMOLOGY.out.gff3.map { meta, gff3 -> gff3 }.first()
    WRITE_MANIFEST(params.outdir, ch_gff3)

    RUN_HOMOLOGY.out.versions.collectFile(
        name: 'software_versions.tsv',
        newLine: true,
        storeDir: "${params.outdir}/pipeline_info"
    )
}
