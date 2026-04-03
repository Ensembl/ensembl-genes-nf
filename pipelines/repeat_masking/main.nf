#!/usr/bin/env nextflow
/*
========================================================================================
    REPEAT_MASKING PIPELINE
========================================================================================
    Replaces the Perl/eHive RepeatMasking subpipeline.

    Steps:
      1. (Optional) RepeatModeler — de novo repeat library construction
      2. RepeatMasker  — interspersed repeat masking (RepBase / custom library)
      3. RED           — rapid de novo repeat detection
      4. TRF           — tandem repeat masking
      5. DUST          — low-complexity masking
      6. Merge all repeat BEDs into a single GFF3
      7. bedtools maskfasta — produce softmasked genome FASTA
      8. Write output_manifest.json for HiveRunNextflow dataflow

    Flat-file I/O: no Ensembl core DB dependency.

    Run locally (stub, no containers):
      nextflow run . -profile local -stub \
        --genome_fasta genome.fa --outdir results/

    Run on Slurm:
      nextflow run . -profile slurm \
        --genome_fasta genome.fa --repbase_library RepBase.h5 --outdir results/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { SAMTOOLS_FAIDX } from './modules/samtools_faidx.nf'
include { WRITE_MANIFEST } from './modules/write_manifest.nf'
include { BUILD_LIBRARY  } from './subworkflows/build_library.nf'
include { MASK_REPEATS   } from './subworkflows/mask_repeats.nf'

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

def validate_params() {
    def errors = []
    if (!params.genome_fasta) errors << "  --genome_fasta is required"
    if (!params.outdir)       errors << "  --outdir is required"
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
    // REFERENCE GENOME
    //
    def genome_meta = [id: file(params.genome_fasta).baseName]
    ch_genome = Channel.of([genome_meta, file(params.genome_fasta, checkIfExists: true)])

    //
    // INDEX GENOME
    //
    SAMTOOLS_FAIDX(ch_genome, [[],[]])

    //
    // SUBWORKFLOW: Build repeat library (optional RepeatModeler + custom)
    //
    BUILD_LIBRARY(
        ch_genome,
        params.skip_repeatmodeler,
        params.custom_library
    )

    //
    // CHUNK GENOME for parallelisation
    // chunk_fasta.py is invoked inside RepeatMasker/RED/TRF/DUST processes
    // via work-dir isolation. For large genomes, replace with a CHUNK_FASTA
    // process that emits one file per chunk.
    //
    ch_chunks = ch_genome
        .map { meta, fasta -> [[meta + [genome_id: meta.id]], fasta] }

    //
    // SUBWORKFLOW: Mask repeats (RepeatMasker + RED + TRF + DUST → merge → mask)
    //
    MASK_REPEATS(
        ch_chunks,
        ch_genome,
        BUILD_LIBRARY.out.library,
        params.skip_red,
        params.skip_trf,
        params.skip_dust
    )

    //
    // WRITE OUTPUT MANIFEST — required by HiveRunNextflow bridge
    //
    WRITE_MANIFEST(
        params.outdir,
        MASK_REPEATS.out.softmasked_fasta.map { meta, fa  -> fa  }.first(),
        MASK_REPEATS.out.repeat_gff3.map      { meta, gff -> gff }.first()
    )

    //
    // SOFTWARE VERSIONS
    //
    ch_versions = Channel.empty()
        .mix(MASK_REPEATS.out.versions)
        .mix(BUILD_LIBRARY.out.versions)
        .mix(SAMTOOLS_FAIDX.out.versions)

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
