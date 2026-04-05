#!/usr/bin/env nextflow
/*
========================================================================================
    UTR_ADDITION PIPELINE
========================================================================================
    Add UTR regions to protein-coding gene models by matching them against
    full-length donor transcripts from RNA-seq, long-read or cDNA alignment pipelines.

    Algorithm overview:
      For each acceptor transcript (from the consolidated geneset — has CDS but
      incomplete or no UTRs):
        1. Extract the internal CDS splice junctions.
        2. Search donor GFF3 files in priority order (first file = highest priority).
           A donor transcript is a match if its exon structure CONTAINS the acceptor
           CDS:  all internal CDS splice sites must appear as exon boundaries in the
           donor; the first/last CDS exon boundaries (which face the UTR) may differ.
        3. Take the first match found (highest-priority donor file wins).
        4. Rebuild the transcript's exon list by replacing the CDS-flanking ends with
           the donor's UTR exons, while keeping the CDS exons exactly as they are.
        5. Clip any UTR that exceeds max_5prime_utr / max_3prime_utr bp.
        6. Discard UTR exons smaller than min_utr_exon_size bp.
      Modified transcripts carry a utr_source= attribute naming the donor file.

    Donor priority:
      Pass donor files in order of decreasing confidence, e.g.:
        long_read.gff3,best_targeted.gff3,rnaseq.gff3

    Inputs:
      --consolidated_gff3   Coding gene models to enrich with UTRs (GFF3)
      --donor_gff3_files    Comma-separated donor GFF3 files (priority order)
      --outdir              Output directory
      --max_5prime_utr      Maximum 5' UTR extension in bp (default: 5000)
      --max_3prime_utr      Maximum 3' UTR extension in bp (default: 10000)
      --min_utr_exon_size   Minimum UTR exon size in bp to retain (default: 30)

    Run locally (stub):
      nextflow run . -profile local -stub \\
        --consolidated_gff3  consolidated.gff3 \\
        --donor_gff3_files   long_read.gff3,rnaseq.gff3 \\
        --outdir             results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { ADD_UTRS       } from './modules/add_utrs.nf'
include { WRITE_MANIFEST } from './modules/write_manifest.nf'

def validate_params() {
    def errors = []
    if (!params.consolidated_gff3)
        errors << "  --consolidated_gff3 is required"
    if (!params.donor_gff3_files)
        errors << "  --donor_gff3_files is required (comma-separated GFF3 paths, priority order)"
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
    // Acceptor gene models (single consolidated GFF3)
    //
    ch_consolidated = file(params.consolidated_gff3, checkIfExists: true)

    //
    // Donor transcripts — collected into a single list preserving input order
    // (order = priority: index 0 is highest priority)
    //
    ch_donors = Channel
        .fromPath(
            params.donor_gff3_files.split(',').collect { it.trim() },
            checkIfExists: true
        )
        .collect()

    //
    // ADD UTRs
    //
    ADD_UTRS(
        ch_consolidated,
        ch_donors,
        params.max_5prime_utr,
        params.max_3prime_utr,
        params.min_utr_exon_size
    )

    //
    // MANIFEST
    //
    WRITE_MANIFEST(params.outdir, ADD_UTRS.out.gff3)

    Channel.empty()
        .mix(ADD_UTRS.out.versions)
        .collectFile(
            name: 'software_versions.tsv',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
