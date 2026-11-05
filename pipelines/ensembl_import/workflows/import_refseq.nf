include { IMPORT_REFSEQ_TO_CORE } from '../subworkflows/import_refseq_to_core.nf'
include { IMPORT_REPEATMASKER } from '../subworkflows/import_repeatmasker.nf'
include { GET_METADATA_CORE } from '../subworkflows/get_metadata_core.nf'
include { PREPARE_STATS_INPUT } from '../modules/prepare_stats_input.nf'
include { COMBINE_STATS_INPUT } from '../modules/combine_stats_input.nf'
include { REGISTRY_TRACKER_START } from '../modules/registry_tracker_start.nf'
include { REGISTRY_TRACKER_COMPLETE } from '../modules/registry_tracker_complete.nf'

workflow IMPORT_REFSEQ {

    take:
        // channel: val(meta)
        //   meta: [id: RefSeq assembly accession, species: binomial species name]
        samples_ch

        // value channel: tuple val(db_host), val(db_port), val(db_user),
        // val(db_password), val(db_read_user)
        db_config_ch

    main:

        REGISTRY_TRACKER_START(samples_ch)

        samples_with_registry_status_ch = REGISTRY_TRACKER_START.out.started
            .map { meta, import_status_file ->
                meta + [registry_import_status_id: import_status_file.text.trim()]
            }

        IMPORT_REFSEQ_TO_CORE(samples_with_registry_status_ch, db_config_ch)

        IMPORT_REPEATMASKER(
            IMPORT_REFSEQ_TO_CORE.out.loaded_refseq,
            IMPORT_REFSEQ_TO_CORE.out.assembly_report,
            db_config_ch
        )

        GET_METADATA_CORE(
            IMPORT_REFSEQ_TO_CORE.out.loaded_refseq,
            db_config_ch
        )

        core_genome_ready = IMPORT_REFSEQ_TO_CORE.out.genome
            .map { meta, genome -> tuple(meta.id, meta, genome) }

        taxonomy_ready = GET_METADATA_CORE.out.taxonomy_loaded
            .map { meta, taxonomy_marker -> tuple(meta.id, taxonomy_marker) }

        repeatmasker_ready = IMPORT_REPEATMASKER.out.completed
            .map { meta, repeatmasker_marker -> tuple(meta.id, repeatmasker_marker) }

        import_ready = core_genome_ready
            .join(taxonomy_ready)
            .join(repeatmasker_ready)
            .map { _id, meta, genome, _taxonomy_marker, _repeatmasker_marker -> tuple(meta, genome) }

        REGISTRY_TRACKER_COMPLETE(
            import_ready.map { meta, _genome -> meta }
        )

        PREPARE_STATS_INPUT(
            import_ready
        )
        combined_stats_input = COMBINE_STATS_INPUT(PREPARE_STATS_INPUT.out.csv.collect())

        versions_ch = IMPORT_REFSEQ_TO_CORE.out.versions
            .mix(IMPORT_REPEATMASKER.out.versions)
            .mix(GET_METADATA_CORE.out.versions)
            .mix(PREPARE_STATS_INPUT.out.versions)
            .mix(COMBINE_STATS_INPUT.out.versions)
            .mix(REGISTRY_TRACKER_START.out.versions)
            .mix(REGISTRY_TRACKER_COMPLETE.out.versions)

    emit:
        loaded_refseq = IMPORT_REFSEQ_TO_CORE.out.loaded_refseq
        stats_input = combined_stats_input.csv
        taxonomy_loaded = GET_METADATA_CORE.out.taxonomy_loaded
        versions = versions_ch
    }
