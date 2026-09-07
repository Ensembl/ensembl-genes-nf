include { RESOLVE_SPECIES_INPUTS } from '../../modules/local/resolve_species_inputs.nf'


workflow RESOLVE_AND_ROUTE {
    main:
    def has_samplesheet = params.samplesheet != null
    def has_db_name = (
        params.db_name != null &&
        params.db_name.toString().trim()
    )

    if (has_samplesheet && has_db_name) {
        error "Provide either --samplesheet or --db_name, not both"
    }

    if (!has_samplesheet && !has_db_name) {
        error "Provide either --samplesheet <file.csv> or --db_name <database>"
    }

    def default_mode = (
        params.mode != null && params.mode.toString().trim()
        ? params.mode.toString().trim()
        : 'auto'
    )

    if (!(default_mode in ['auto', 'map', 'reassign'])) {
        error "Invalid --mode '${default_mode}'; expected auto, map, or reassign"
    }

    def direct_overrides = [
        params.target_fasta,
        params.target_gff,
        params.ref_fasta,
        params.ref_gff,
        params.mapping_session_id,
    ]

    if (
        has_samplesheet &&
        direct_overrides.any { value ->
            value != null && value.toString().trim()
        }
    ) {
        error "Direct input overrides require --db_name and cannot be used with --samplesheet"
    }

    if (has_samplesheet) {
        species_ch = Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def db_name = (
                    row.get('db_name') != null
                    ? row.get('db_name').toString().trim()
                    : ''
                )

                if (!db_name) {
                    error "Every samplesheet row must contain db_name"
                }

                def row_mode = (
                    row.get('mode') != null
                    ? row.get('mode').toString().trim()
                    : ''
                )
                def requested_mode = row_mode ?: default_mode

                if (!(requested_mode in ['auto', 'map', 'reassign'])) {
                    error(
                        "Invalid mode '${requested_mode}' for ${db_name}; " +
                        "expected auto, map, or reassign"
                    )
                }

                tuple(
                    db_name,
                    requested_mode,
                    row.get('target_fasta') ?: '',
                    row.get('target_gff') ?: '',
                    row.get('mapping_session_id') ?: '',
                    row.get('ref_fasta') ?: '',
                    row.get('ref_gff') ?: ''
                )
            }
    } else {
        species_ch = Channel.of(
            tuple(
                params.db_name.toString().trim(),
                default_mode,
                params.target_fasta ?: '',
                params.target_gff ?: '',
                params.mapping_session_id ?: '',
                params.ref_fasta ?: '',
                params.ref_gff ?: ''
            )
        )
    }

    RESOLVE_SPECIES_INPUTS(species_ch)

    resolved_routes = RESOLVE_SPECIES_INPUTS.out.inputs_json
        .map { json_file ->
            def data = new groovy.json.JsonSlurper().parse(json_file)

            if (data.effective_mode == 'no_action') {
                log.info(
                    "NO ACTION REQUIRED: ${data.db_name} stable IDs " +
                    "agree with the registry allocation."
                )
            }

            data
        }
        .branch {
            mapping:
                it.effective_mode == 'map'

            reassignment:
                it.effective_mode == 'reassign'
        }

    mapping_ch = resolved_routes.mapping.map { data ->
        tuple(
            data.db_name,
            file(data.ref_fasta),
            file(data.ref_gff),
            file(data.target_fasta),
            file(data.target_gff),
            data.mapping_session_id,
            data.gene_range,
            data.transcript_range,
            data.translation_range
        )
    }

    reassignment_ch = resolved_routes.reassignment.map { data ->
        tuple(
            data.db_name,
            data.gene_range,
            data.transcript_range,
            data.translation_range
        )
    }

    emit:
    mapping = mapping_ch
    reassignment = reassignment_ch
}
