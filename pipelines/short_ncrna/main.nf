#!/usr/bin/env nextflow
/*
========================================================================================
    SHORT_NCRNA PIPELINE
========================================================================================
    Short non-coding RNA annotation via Rfam (cmsearch) + miRBase (BLASTN).
    Replaces the Perl/eHive ShortncRNA subpipeline.

    Steps:
      1. Split genome into 1 Mb chunks for parallel search
      2. cmsearch with Rfam covariance models (--cut_ga for per-family thresholds)
      3. (Optional) BLASTN of miRBase sequences against genome
      4. Filter by e-value/score, assign biotypes from CM names
      5. Merge all chunk GFF3 → final ncRNA annotation
      6. Write output_manifest.json for HiveRunNextflow dataflow

    Input genome: UNMASKED (unlike repeat_masking; cmsearch performs better unmasked)
    Biotypes produced: miRNA, snRNA, snoRNA, scaRNA, rRNA, tRNA, ribozyme,
                       vault_RNA, Y_RNA, SRP_RNA, misc_RNA

    Run locally (stub):
      nextflow run . -profile local -stub \
        --genome_fasta genome.fa --rfam_cm Rfam.cm --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --genome_fasta genome.fa --rfam_cm Rfam.cm \
        [--mirna_fasta all_mirnas.fa] --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { WRITE_MANIFEST } from './modules/write_manifest.nf'
include { SEARCH_NCRNA   } from './subworkflows/search_ncrna.nf'

def validate_params() {
    def errors = []
    if (!params.genome_fasta) errors << "  --genome_fasta is required (UNMASKED genome)"
    if (!params.rfam_cm)      errors << "  --rfam_cm is required (Rfam covariance models)"
    if (!params.outdir)       errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    //
    // CHUNK GENOME into 1 Mb slices for parallel cmsearch
    //
    ch_genome_chunks = Channel
        .fromPath(params.genome_fasta, checkIfExists: true)
        .splitFasta(by: params.chunk_size, file: true)
        .map { fa ->
            def chunk_id = fa.name.replaceAll(/\.fasta$|\.fa$/, '')
            [[id: chunk_id], fa]
        }

    //
    // SUBWORKFLOW: cmsearch + BLASTN per chunk
    //
    SEARCH_NCRNA(
        ch_genome_chunks,
        file(params.rfam_cm, checkIfExists: true),
        params.mirna_fasta ?: null,
        params.mirna_blast_db ? [[id: 'genome_db'], file(params.mirna_blast_db, type: 'dir')] : null
    )

    //
    // MERGE all chunk GFF3 files into a single output
    //
    ch_all_gff3 = SEARCH_NCRNA.out.gff3
        .map { meta, gff3 -> gff3 }
        .collect()
        .map { files -> files }

    // Collect to single file for manifest
    ch_merged = ch_all_gff3.map { files -> files[0] }.first()

    WRITE_MANIFEST(params.outdir, ch_merged)

    SEARCH_NCRNA.out.versions.collectFile(
        name: 'software_versions.tsv',
        newLine: true,
        storeDir: "${params.outdir}/pipeline_info"
    )
}
