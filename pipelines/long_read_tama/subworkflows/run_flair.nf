nextflow.enable.dsl = 2

include { FLAIR_JUNCTIONS; FLAIR_TRANSCRIPTOME; FLAIR_COMBINE } from '../modules/flair.nf'
include { GTF_TO_BED12 } from '../modules/gtf_to_bed12.nf'

workflow RUN_FLAIR {
    take:
    aligned_bams
    reference

    main:
    // FLAIR transcriptome is run on complete accession BAMs. Its native
    // combine step is retained for the cohort result; no other backend's
    // merge tool is used.
    FLAIR_JUNCTIONS(aligned_bams)
    flair_inputs = aligned_bams
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(FLAIR_JUNCTIONS.out.bed.map { meta, junctions -> tuple(meta.id, junctions) })
        .map { _id, meta, bam, bai, junctions -> tuple(meta, bam, bai, junctions) }
    FLAIR_TRANSCRIPTOME(flair_inputs, reference)
    accession_products = FLAIR_TRANSCRIPTOME.out.products
    cohort_input = accession_products
        .map { _meta, archive -> archive }
        .collect()
        .map { archives -> tuple([id: params.cohort_id, scope: 'cohort'], archives) }
    FLAIR_COMBINE(cohort_input)
    accession_gtf = FLAIR_TRANSCRIPTOME.out.gtf.map { meta, gtf -> tuple(meta, 'flair', 'accession', gtf) }
    cohort_gtf = FLAIR_COMBINE.out.gtf.map { meta, gtf -> tuple(meta, 'flair', 'cohort', gtf) }
    GTF_TO_BED12(accession_gtf.mix(cohort_gtf))

    emit:
    products = accession_products.mix(FLAIR_COMBINE.out.products)
    bed = GTF_TO_BED12.out.bed.filter { _meta, _backend, scope, _bed -> scope == 'cohort' }
        .map { meta, _backend, scope, bed -> tuple(meta, 'flair', scope, bed) }
    accession_bed = GTF_TO_BED12.out.bed.filter { _meta, _backend, scope, _bed -> scope == 'accession' }
        .map { meta, _backend, scope, bed -> tuple(meta, 'flair', scope, bed) }
    reports = accession_products.mix(FLAIR_COMBINE.out.products)
    versions = FLAIR_TRANSCRIPTOME.out.versions
        .mix(FLAIR_JUNCTIONS.out.versions)
        .mix(FLAIR_COMBINE.out.versions)
        .mix(GTF_TO_BED12.out.versions)
}
