include { INPUT_PREPARATION } from '../subworkflows/input_preparation.nf'
include { SUBSTRATE_TOPOLOGY_WORKFLOW } from '../subworkflows/substrate_topology.nf'
include { TRANSLATION_RECONCILIATION } from '../subworkflows/translation_reconciliation.nf'
include { INSTANCE_CONTEXT } from '../subworkflows/instance_context.nf'
include { CONSTRAINT_AXES } from '../subworkflows/constraint_axes.nf'
include { PEPTIDE_PREPARATION } from '../subworkflows/peptide_preparation.nf'
include { PEPTIDE_AXES } from '../subworkflows/peptide_axes.nf'
include { OBSERVED_EVIDENCE } from '../subworkflows/observed_evidence.nf'
include { MERGE_TYPED_AXES } from '../subworkflows/merge_axes.nf'
include { PRODUCT_AND_ADJUDICATION } from '../subworkflows/product_and_adjudication.nf'

workflow TRANSLON_CHARACTERISATION {
    take:
    intervals
    gencode_gff3
    translation_verdicts
    proteome_fasta
    genome_fasta

    main:
    INPUT_PREPARATION(intervals, gencode_gff3, translation_verdicts, proteome_fasta, genome_fasta)
    SUBSTRATE_TOPOLOGY_WORKFLOW(INPUT_PREPARATION.out.intervals, INPUT_PREPARATION.out.annotation, INPUT_PREPARATION.out.genome)
    TRANSLATION_RECONCILIATION(SUBSTRATE_TOPOLOGY_WORKFLOW.out.instances, INPUT_PREPARATION.out.verdicts)

    // Independent per-instance work: each subworkflow preserves the same meta tuple.
    INSTANCE_CONTEXT(TRANSLATION_RECONCILIATION.out.instances)
    CONSTRAINT_AXES(TRANSLATION_RECONCILIATION.out.instances, INPUT_PREPARATION.out.genome)
    PEPTIDE_PREPARATION(TRANSLATION_RECONCILIATION.out.instances)
    PEPTIDE_AXES(PEPTIDE_PREPARATION.out.peptides, INPUT_PREPARATION.out.proteome, TRANSLATION_RECONCILIATION.out.instances)
    OBSERVED_EVIDENCE(PEPTIDE_PREPARATION.out.peptides)

    all_axes = INSTANCE_CONTEXT.out.axes
        .concat(CONSTRAINT_AXES.out.axes)
        .concat(PEPTIDE_AXES.out.axes)
        .concat(OBSERVED_EVIDENCE.out.axes)
    MERGE_TYPED_AXES(all_axes)
    PRODUCT_AND_ADJUDICATION(MERGE_TYPED_AXES.out.instances)

    emit:
    unhosted = SUBSTRATE_TOPOLOGY_WORKFLOW.out.unhosted
    claims = PRODUCT_AND_ADJUDICATION.out.claims
    clusters = PRODUCT_AND_ADJUDICATION.out.clusters
    manifest = PRODUCT_AND_ADJUDICATION.out.manifest
    versions = SUBSTRATE_TOPOLOGY_WORKFLOW.out.versions
        .concat(INPUT_PREPARATION.out.versions)
        .concat(TRANSLATION_RECONCILIATION.out.versions)
        .concat(INSTANCE_CONTEXT.out.versions)
        .concat(CONSTRAINT_AXES.out.versions)
        .concat(PEPTIDE_AXES.out.versions)
        .concat(PEPTIDE_PREPARATION.out.versions)
        .concat(OBSERVED_EVIDENCE.out.versions)
        .concat(MERGE_TYPED_AXES.out.versions)
        .concat(PRODUCT_AND_ADJUDICATION.out.versions)
}
