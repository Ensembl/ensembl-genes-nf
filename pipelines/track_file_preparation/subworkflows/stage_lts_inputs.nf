include { STAGE_LTS_FILES } from '../modules/stage_lts_files.nf'

workflow STAGE_LTS_INPUTS {
    take:
    entities

    main:
    branched = entities.branch { _meta, _bam, _gtf, _sj, _bai ->
        needs_staging: [_bam, _gtf, _sj, _bai].any { input -> input && input.toString().startsWith('/lts/') }
        direct: true
    }

    staging_inputs = branched.needs_staging.map { meta, bam, gtf, sj, bai ->
        def fields = [bam: bam, gtf: gtf, sj_out_tab: sj, bai: bai]
        def lts_fields = fields.findAll { _name, input -> input && input.toString().startsWith('/lts/') }
        def stage_kinds = lts_fields.keySet().toList()
        def sources = lts_fields.values().toList()
        def stage_names = stage_kinds.withIndex().collect { kind, index ->
            def basename = sources[index].getFileName().toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
            "stage_${String.format('%02d', index + 1)}__${meta.gca_accession}__${meta.safe_id}__${kind}__${basename}"
        }
        def source_fields = fields.collectEntries { _name, input -> [(_name): input?.toString()] }
        tuple(meta, source_fields, stage_kinds, stage_names, sources)
    }
    STAGE_LTS_FILES(staging_inputs)

    localized_entities = STAGE_LTS_FILES.out.staged.map { meta, source_fields, stage_kinds, staged_files ->
        def localized = source_fields.collectEntries { name, input -> [(name): input ? file(input) : null] }
        stage_kinds.eachWithIndex { kind, index ->
            def field_name = kind == 'sj_out_tab' ? 'sj_out_tab' : kind
            localized[field_name] = staged_files[index]
        }
        tuple(meta, localized.bam, localized.gtf, localized.sj_out_tab, localized.bai)
    }

    direct_entities = branched.direct
    localized_entities = localized_entities.mix(direct_entities)

    emit:
    entities = localized_entities
    versions = STAGE_LTS_FILES.out.versions
}
