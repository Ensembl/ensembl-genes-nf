nextflow.enable.dsl = 2

include { GTF_TO_BED12 } from '../modules/gtf_to_bed12.nf'
include { GTF_TO_BED12 as GTF_TO_BED12_ACCESSION } from '../modules/gtf_to_bed12.nf'
include { AUDIT_NATIVE_MODELS; AUDIT_NATIVE_MODELS as AUDIT_NATIVE_COHORT_MODELS } from '../modules/audit_native_models.nf'
include { STRINGTIE2_MERGE; STRINGTIE2_MERGE as STRINGTIE2_COHORT_MERGE } from '../modules/stringtie2_merge.nf'
include { STRINGTIE3_MERGE; STRINGTIE3_MERGE as STRINGTIE3_COHORT_MERGE } from '../modules/stringtie3_merge.nf'
include { TMERGE_NATIVE_MERGE; TMERGE_NATIVE_MERGE as TMERGE_NATIVE_COHORT_MERGE } from '../modules/tmerge_native_merge.nf'
include { TAMA_MERGE_ACCESSION } from '../modules/tama_merge_accession.nf'
include { TMERGE_MERGE_ACCESSION } from '../modules/tmerge_merge_accession.nf'
include { TAMA_MERGE } from '../modules/tama_merge.nf'
include { TMERGE } from '../modules/tmerge.nf'

workflow MERGE_LONG_READ_MODELS {
    take:
    native_models

    main:
    // native_models: tuple val(meta), val(backend), val(shard),
    //                path(native_model), val(native_format)
    native_mode = params.backend_merge_mode == 'native'
    native_cohort_models = channel.empty()
    native_reports = channel.empty()
    native_versions = channel.empty()
    legacy_beds = channel.empty()
    accession_beds = channel.empty()
    legacy_reports = channel.empty()
    legacy_versions = channel.empty()

    if (native_mode) {
        accession_groups = native_models
            .map { meta, backend, shard, model, native_format -> tuple("${backend}@@${meta.id}", meta, backend, shard, model, native_format) }
            .groupTuple()

        tama_accession_inputs = accession_groups.filter { _key, _metas, backends, _shards, _models, _formats -> backends[0] == 'tama' }
            .map { _key, metas, backends, _shards, models, _formats -> tuple(metas[0], backends[0], metas[0].id, models.sort { left, right -> left.name <=> right.name }) }
        stringtie2_accession_inputs = accession_groups.filter { _key, _metas, backends, _shards, _models, _formats -> backends[0] == 'stringtie2' }
            .map { _key, metas, backends, _shards, models, _formats -> tuple(metas[0], backends[0], metas[0].id, models.sort { left, right -> left.name <=> right.name }) }
        stringtie3_accession_inputs = accession_groups.filter { _key, _metas, backends, _shards, _models, _formats -> backends[0] == 'stringtie3' }
            .map { _key, metas, backends, _shards, models, _formats -> tuple(metas[0], backends[0], metas[0].id, models.sort { left, right -> left.name <=> right.name }) }
        tmerge_accession_inputs = accession_groups.filter { _key, _metas, backends, _shards, _models, _formats -> backends[0] == 'tmerge' }
            .map { _key, metas, backends, _shards, models, _formats -> tuple(metas[0], backends[0], metas[0].id, models.sort { left, right -> left.name <=> right.name }) }

        TAMA_MERGE_ACCESSION(tama_accession_inputs)
        STRINGTIE2_MERGE(stringtie2_accession_inputs)
        STRINGTIE3_MERGE(stringtie3_accession_inputs)
        TMERGE_NATIVE_MERGE(tmerge_accession_inputs)

        accession_native_models = TAMA_MERGE_ACCESSION.out.bed
            .map { meta, backend, accession, bed -> tuple(meta, backend, accession, bed, 'bed12') }
            .mix(STRINGTIE2_MERGE.out.gtf.map { meta, backend, accession, gtf -> tuple(meta, backend, accession, gtf, 'gtf') })
            .mix(STRINGTIE3_MERGE.out.gtf.map { meta, backend, accession, gtf -> tuple(meta, backend, accession, gtf, 'gtf') })
            .mix(TMERGE_NATIVE_MERGE.out.gtf.map { meta, backend, accession, gtf -> tuple(meta, backend, accession, gtf, 'gtf') })

        AUDIT_NATIVE_MODELS(accession_native_models.map { meta, backend, accession, model, native_format -> tuple(meta, backend, accession, model, native_format) })

        cohort_groups = accession_native_models
            .map { meta, backend, accession, model, native_format -> tuple("${backend}@@${params.cohort_id}", meta, backend, accession, model, native_format) }
            .groupTuple()

        tama_cohort_inputs = cohort_groups.filter { _key, _metas, backends, _accessions, _models, _formats -> backends[0] == 'tama' }
            .map { _key, metas, backends, _accessions, models, _formats -> tuple(metas[0] + [id: params.cohort_id], backends[0], params.cohort_id, models.sort { left, right -> left.name <=> right.name }) }
        stringtie2_cohort_inputs = cohort_groups.filter { _key, _metas, backends, _accessions, _models, _formats -> backends[0] == 'stringtie2' }
            .map { _key, metas, backends, _accessions, models, _formats -> tuple(metas[0] + [id: params.cohort_id], backends[0], params.cohort_id, models.sort { left, right -> left.name <=> right.name }) }
        stringtie3_cohort_inputs = cohort_groups.filter { _key, _metas, backends, _accessions, _models, _formats -> backends[0] == 'stringtie3' }
            .map { _key, metas, backends, _accessions, models, _formats -> tuple(metas[0] + [id: params.cohort_id], backends[0], params.cohort_id, models.sort { left, right -> left.name <=> right.name }) }
        tmerge_cohort_inputs = cohort_groups.filter { _key, _metas, backends, _accessions, _models, _formats -> backends[0] == 'tmerge' }
            .map { _key, metas, backends, _accessions, models, _formats -> tuple(metas[0] + [id: params.cohort_id], backends[0], params.cohort_id, models.sort { left, right -> left.name <=> right.name }) }

        TAMA_MERGE(tama_cohort_inputs)
        STRINGTIE2_COHORT_MERGE(stringtie2_cohort_inputs)
        STRINGTIE3_COHORT_MERGE(stringtie3_cohort_inputs)
        TMERGE_NATIVE_COHORT_MERGE(tmerge_cohort_inputs)

        native_cohort_models = TAMA_MERGE.out.bed
            .map { meta, backend, cohort, bed -> tuple(meta, backend, cohort, bed) }
            .mix(STRINGTIE2_COHORT_MERGE.out.gtf.map { meta, backend, cohort, gtf -> tuple(meta, backend, cohort, gtf) })
            .mix(STRINGTIE3_COHORT_MERGE.out.gtf.map { meta, backend, cohort, gtf -> tuple(meta, backend, cohort, gtf) })
            .mix(TMERGE_NATIVE_COHORT_MERGE.out.gtf.map { meta, backend, cohort, gtf -> tuple(meta, backend, cohort, gtf) })

        AUDIT_NATIVE_COHORT_MODELS(native_cohort_models.map { meta, backend, cohort, model -> tuple(meta, backend, cohort, model, model.name.endsWith('.gtf') ? 'gtf' : 'bed12') })

        native_gtf_inputs = native_cohort_models.filter { _meta, _backend, _cohort, model -> model.name.endsWith('.gtf') }
            .map { meta, backend, cohort, gtf -> tuple(meta, backend, cohort, gtf) }
        GTF_TO_BED12(native_gtf_inputs)
        native_cohort_models = native_cohort_models.filter { _meta, _backend, _cohort, model -> model.name.endsWith('.bed') }
            .mix(GTF_TO_BED12.out.bed.map { meta, backend, cohort, bed -> tuple(meta, backend, cohort, bed) })

        native_accession_gtf_inputs = accession_native_models.filter { _meta, _backend, _accession, _model, native_format -> native_format == 'gtf' }
            .map { meta, backend, accession, gtf, _format -> tuple(meta, backend, accession, gtf) }
        GTF_TO_BED12_ACCESSION(native_accession_gtf_inputs)
        accession_beds = accession_native_models.filter { _meta, _backend, _accession, _model, native_format -> native_format == 'bed12' }
            .map { meta, backend, accession, bed, _format -> tuple(meta, backend, accession, bed) }
            .mix(GTF_TO_BED12_ACCESSION.out.bed)

        native_reports = TAMA_MERGE_ACCESSION.out.filelist.mix(TAMA_MERGE_ACCESSION.out.filelist_checksum)
            .mix(AUDIT_NATIVE_MODELS.out.manifest).mix(AUDIT_NATIVE_MODELS.out.stats).mix(AUDIT_NATIVE_MODELS.out.checksum)
            .mix(STRINGTIE2_MERGE.out.filelist).mix(STRINGTIE2_MERGE.out.filelist_checksum)
            .mix(STRINGTIE3_MERGE.out.filelist).mix(STRINGTIE3_MERGE.out.filelist_checksum)
            .mix(TMERGE_NATIVE_MERGE.out.filelist).mix(TMERGE_NATIVE_MERGE.out.filelist_checksum)
            .mix(STRINGTIE2_COHORT_MERGE.out.filelist).mix(STRINGTIE2_COHORT_MERGE.out.filelist_checksum)
            .mix(STRINGTIE3_COHORT_MERGE.out.filelist).mix(STRINGTIE3_COHORT_MERGE.out.filelist_checksum)
            .mix(TMERGE_NATIVE_COHORT_MERGE.out.filelist).mix(TMERGE_NATIVE_COHORT_MERGE.out.filelist_checksum)
            .mix(TAMA_MERGE.out.filelist).mix(TAMA_MERGE.out.filelist_checksum)
            .mix(AUDIT_NATIVE_COHORT_MODELS.out.manifest).mix(AUDIT_NATIVE_COHORT_MODELS.out.stats).mix(AUDIT_NATIVE_COHORT_MODELS.out.checksum)
            .mix(GTF_TO_BED12_ACCESSION.out.versions)
        native_versions = TAMA_MERGE_ACCESSION.out.versions
            .mix(STRINGTIE2_MERGE.out.versions).mix(STRINGTIE3_MERGE.out.versions).mix(TMERGE_NATIVE_MERGE.out.versions)
            .mix(TAMA_MERGE.out.versions).mix(STRINGTIE2_COHORT_MERGE.out.versions).mix(STRINGTIE3_COHORT_MERGE.out.versions).mix(TMERGE_NATIVE_COHORT_MERGE.out.versions)
            .mix(GTF_TO_BED12.out.versions).mix(GTF_TO_BED12_ACCESSION.out.versions)
    } else {
        legacy_gtf_inputs = native_models.filter { _meta, _backend, _shard, _model, native_format -> native_format == 'gtf' }
            .map { meta, backend, shard, gtf, _format -> tuple(meta, backend, shard, gtf) }
        GTF_TO_BED12(legacy_gtf_inputs)
        legacy_beds = native_models.filter { _meta, _backend, _shard, _model, native_format -> native_format == 'bed12' }
            .map { meta, backend, shard, bed, _format -> tuple(meta, backend, shard, bed) }
            .mix(GTF_TO_BED12.out.bed)

        accession_inputs = legacy_beds.map { meta, backend, _shard, bed -> tuple("${backend}@@${meta.id}", meta, backend, meta.id, bed) }
            .groupTuple().map { _key, metas, backends, accessions, beds -> tuple(metas[0], backends[0], accessions[0], beds.sort { left, right -> left.name <=> right.name }) }

        if (params.merge_tool == 'tmerge') {
            TMERGE_MERGE_ACCESSION(accession_inputs)
            accession_beds = TMERGE_MERGE_ACCESSION.out.bed
            legacy_reports = TMERGE_MERGE_ACCESSION.out.report
            legacy_versions = TMERGE_MERGE_ACCESSION.out.versions
        } else {
            TAMA_MERGE_ACCESSION(accession_inputs)
            accession_beds = TAMA_MERGE_ACCESSION.out.bed
            legacy_reports = TAMA_MERGE_ACCESSION.out.report
            legacy_versions = TAMA_MERGE_ACCESSION.out.versions
        }

        cohort_inputs = accession_beds.map { meta, backend, _accession, bed -> tuple("${backend}@@${params.cohort_id}", meta, backend, params.cohort_id, bed) }
            .groupTuple().map { _key, metas, backends, cohorts, models -> tuple(metas[0] + [id: params.cohort_id], backends[0], cohorts[0], models.sort { left, right -> left.name <=> right.name }) }
        if (params.merge_tool == 'tmerge') {
            TMERGE(cohort_inputs)
            legacy_beds = TMERGE.out.bed.map { meta, backend, cohort, bed -> tuple(meta, backend, cohort, bed) }
            legacy_reports = legacy_reports.mix(TMERGE.out.merge_report).mix(TMERGE.out.trans_report)
            legacy_versions = legacy_versions.mix(TMERGE.out.versions)
        } else {
            TAMA_MERGE(cohort_inputs)
            legacy_beds = TAMA_MERGE.out.bed.map { meta, backend, cohort, bed -> tuple(meta, backend, cohort, bed) }
            legacy_reports = legacy_reports.mix(TAMA_MERGE.out.gene_report).mix(TAMA_MERGE.out.merge_report).mix(TAMA_MERGE.out.trans_report)
            legacy_versions = legacy_versions.mix(TAMA_MERGE.out.versions)
        }
        legacy_versions = legacy_versions.mix(GTF_TO_BED12.out.versions)
    }

    emit:
    bed = native_mode ? native_cohort_models : legacy_beds
    accession_bed = accession_beds
    reports = native_mode ? native_reports : legacy_reports
    versions = native_mode ? native_versions : legacy_versions
}
