include { TYPED_AXIS } from '../modules/typed_axis.nf'
include { TYPED_AXIS as PHYLOCSF_UNINFORMATIVE } from '../modules/typed_axis.nf'
include { TYPED_AXIS as MSA_UNINFORMATIVE } from '../modules/typed_axis.nf'
include { PREPARE_MSA_INPUTS } from '../modules/prepare_msa_inputs.nf'
include { ORBL_TOOLS } from '../modules/orbl_tools.nf'
include { FANBACK_ORBL } from '../modules/fanback_orbl.nf'
include { ORBL_ALIGNMENT_DOWNLOAD } from '../modules/orbl_alignment_download.nf'
include { PHYLOCSF_SCORE } from '../modules/phylocsf_score.nf'
include { PHYLOCSF_PARSE } from '../modules/phylocsf_parse.nf'
include { FANBACK_PHYLOCSF } from '../modules/fanback_phylocsf.nf'
include { MSA_REFERENCE_SETUP } from '../modules/msa_reference_setup.nf'
include { MSA_EXTRACT } from '../modules/msa_extract.nf'
include { FANBACK_MSA } from '../modules/fanback_msa.nf'

workflow CONSTRAINT_AXES {
    take:
    instances // [meta, reconciled instances]
    genome // [meta, reference FASTA]

    main:
    // Emit exact instance exon blocks now; MSA_EXTRACT joins these to chromosome MAFs.
    PREPARE_MSA_INPUTS(instances)
    if (params.maf_reference_manifest) {
        MSA_REFERENCE_SETUP(file(params.maf_reference_manifest, checkIfExists: true), instances, genome)
        MSA_EXTRACT(instances, MSA_REFERENCE_SETUP.out.maf_dir)
        FANBACK_MSA(instances, MSA_EXTRACT.out.axis)
        msa_axes = FANBACK_MSA.out.axis.map { meta, file -> tuple(meta, 'msa_extract', file) }
        msa_axis_versions = MSA_EXTRACT.out.versions.concat(FANBACK_MSA.out.versions)
        msa_reference_versions = MSA_REFERENCE_SETUP.out.versions.concat(MSA_EXTRACT.out.versions).concat(FANBACK_MSA.out.versions)
    } else if (params.maf_dir) {
        MSA_EXTRACT(instances, channel.value(file(params.maf_dir, checkIfExists: true)))
        FANBACK_MSA(instances, MSA_EXTRACT.out.axis)
        msa_axes = FANBACK_MSA.out.axis.map { meta, file -> tuple(meta, 'msa_extract', file) }
        msa_axis_versions = MSA_EXTRACT.out.versions.concat(FANBACK_MSA.out.versions)
        msa_reference_versions = MSA_EXTRACT.out.versions.concat(FANBACK_MSA.out.versions)
    } else {
        msa_requests = instances.map { meta, file -> tuple(meta, 'msa_extract', file) }
        MSA_UNINFORMATIVE(msa_requests)
        msa_axes = MSA_UNINFORMATIVE.out.axis
        msa_axis_versions = MSA_UNINFORMATIVE.out.versions
        msa_reference_versions = channel.empty()
    }
    ORBL_TOOLS(instances)
    ORBL_ALIGNMENT_DOWNLOAD(instances)
    FANBACK_ORBL(instances, ORBL_TOOLS.out.axis)
    if (params.phylocsf_matched_null) {
        matched_null = channel.fromPath(params.phylocsf_matched_null, checkIfExists: true)
        PHYLOCSF_SCORE(MSA_EXTRACT.out.alignments)
        phylocsf_parse_input = PHYLOCSF_SCORE.out.scored.combine(matched_null)
            .map { meta, identities, raw, null_file -> tuple(meta, identities, raw, null_file) }
        PHYLOCSF_PARSE(phylocsf_parse_input)
        FANBACK_PHYLOCSF(instances, PHYLOCSF_PARSE.out.axis)
        phylocsf_axes = FANBACK_PHYLOCSF.out.axis.map { meta, file -> tuple(meta, 'phylocsf', file) }
        phylocsf_versions = PHYLOCSF_SCORE.out.versions.concat(PHYLOCSF_PARSE.out.versions).concat(FANBACK_PHYLOCSF.out.versions)
        constraint_names = channel.of('pop_constraint', 'gpn_score')
    } else {
        phylocsf_requests = instances.map { meta, file -> tuple(meta, 'phylocsf', file) }
        PHYLOCSF_UNINFORMATIVE(phylocsf_requests)
        phylocsf_axes = PHYLOCSF_UNINFORMATIVE.out.axis
        phylocsf_versions = PHYLOCSF_UNINFORMATIVE.out.versions
        constraint_names = channel.of('pop_constraint', 'gpn_score')
    }
    requests = instances.combine(constraint_names)
        .map { meta, file, axis -> tuple(meta, axis, file) }
    TYPED_AXIS(requests)

    emit:
    axes = TYPED_AXIS.out.axis.concat(msa_axes).concat(phylocsf_axes).concat(FANBACK_ORBL.out.axis.map { meta, file -> tuple(meta, 'orbl', file) })
    msa_manifest = PREPARE_MSA_INPUTS.out.manifest
    downloaded_alignments = ORBL_ALIGNMENT_DOWNLOAD.out.alignments
    versions = TYPED_AXIS.out.versions.concat(msa_axis_versions).concat(phylocsf_versions).concat(msa_reference_versions).concat(PREPARE_MSA_INPUTS.out.versions).concat(ORBL_TOOLS.out.versions).concat(ORBL_ALIGNMENT_DOWNLOAD.out.versions).concat(FANBACK_ORBL.out.versions)
}
