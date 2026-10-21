#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'
include { PREFLIGHT_TRACK_INPUTS } from './subworkflows/preflight.nf'
include { STAGE_LTS_INPUTS } from './subworkflows/stage_lts_inputs.nf'
include { PREPARE_COVERAGE_TRACKS } from './subworkflows/coverage.nf'
include { INDEX_BAM } from './modules/index_bam.nf'
include { PREPARE_GENE_MODEL_TRACKS } from './subworkflows/gene_models.nf'
include { PREPARE_SPLICE_TRACKS } from './subworkflows/splice.nf'
include { ASSEMBLE_TRACK_SETS } from './subworkflows/reporting.nf'

def safe_entity_id(String value) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def is_lts_location(String raw, base_dir) {
    if (!raw?.trim()) return false
    def path = raw.startsWith('/') ? raw : new File(base_dir.toString(), raw).getAbsolutePath()
    path == '/lts' || path.startsWith('/lts/')
}

def parse_track_selection(String value) {
    def selected = value?.trim() ? value.split(',')*.trim() : ['all']
    if (selected == ['all'])
        return ['coverage', 'gene_model', 'splice_junction']
    def allowed = ['coverage', 'gene_model', 'splice_junction']
    if (selected.any { !(it in allowed) })
        error "Unsupported track_types '${value}'; choose coverage, gene_model, splice_junction, or all"
    selected.unique()
}

def resolve_manifest_path(String raw, base_dir, boolean must_exist, int row_number) {
    if (!raw?.trim()) {
        if (must_exist) error "Manifest row ${row_number} requires a path"
        return null
    }
    def path = file(raw.startsWith('/') ? raw : "${base_dir}/${raw}")
    if (must_exist && !path.exists()) error "Manifest row ${row_number} path does not exist: ${path}"
    path.exists() ? path : null
}

def resolve_current_bam_index(bam, supplied_index, String entity_id, boolean warn_on_stale = true) {
    if (!bam) return null
    def candidates = []
    if (supplied_index) candidates << supplied_index
    [file("${bam}.bai"), file("${bam}.csi")].each { candidate ->
        if (candidate.exists() && !candidates.contains(candidate)) candidates << candidate
    }
    if (!candidates) return null

    def bam_mtime = bam.toFile().lastModified()
    def current = candidates.findAll { it.toFile().lastModified() >= bam_mtime }
    if (warn_on_stale && current.size() < candidates.size())
        log.warn "Ignoring stale BAM index for ${entity_id}: BAM file is newer than its index"
    current.sort { left, right ->
        def left_bai = left.toString().endsWith('.bai') ? 0 : 1
        def right_bai = right.toString().endsWith('.bai') ? 0 : 1
        left_bai <=> right_bai
    }
    current ? current[0] : null
}

def derive_junctions_enabled() {
    params.derive_junctions == true || params.derive_junctions?.toString()?.toBoolean()
}

def manifest_entities(manifest, base_dir, override) {
    def manifest_file = new File(manifest.toString()).getAbsoluteFile()
    def manifest_base = base_dir ?: manifest_file.getParentFile()
    def lines = manifest_file.readLines()
    if (!lines)
        error "Input manifest is empty: ${manifest}"
    def header = lines[0].split('\t', -1)
    def required = ['entity_id', 'entity_type', 'gca_accession', 'assembly_release', 'bam', 'gtf', 'track_types']
    def missing = required.findAll { !(it in header) }
    if (missing)
        error "Manifest is missing required columns: ${missing.join(', ')}"
    def indexes = header.toList().withIndex().collectEntries { pair -> [(pair[0]): pair[1]] }
    def seen = [] as Set
    def safe_seen = [] as Set
    def records = []
    lines.drop(1).eachWithIndex { line, row_number ->
        if (!line.trim()) return
        def values = line.split('\t', -1)
        if (values.size() != header.size())
            error "Manifest row ${row_number + 2} has ${values.size()} fields; expected ${header.size()}"
        def row = header.collectEntries { key -> [(key): values[indexes[key]]] }
        if (!(row.entity_type in ['run', 'sample', 'merged']))
            error "Manifest row ${row_number + 2} has unsupported entity_type '${row.entity_type}'"
        if (!(row.gca_accession ==~ /GCA_[0-9]+\.[0-9]+/))
            error "Manifest row ${row_number + 2} has malformed gca_accession '${row.gca_accession}'"
        def key = "${row.gca_accession}\t${row.entity_id}"
        if (!seen.add(key)) error "Duplicate entity_id plus gca_accession: ${key}"
        def safe = safe_entity_id(row.entity_id)
        if (!safe_seen.add("${row.gca_accession}\t${safe}")) error "Entity ID sanitization collision for ${key}"
        def contains_lts_path = ['bam', 'gtf', 'sj_out_tab', 'bai'].any { field -> is_lts_location(row[field], manifest_base) }
        def selected = override ? parse_track_selection(override) : parse_track_selection(row.track_types)
        def bam = resolve_manifest_path(row.bam, manifest_base, 'coverage' in selected || derive_junctions_enabled(), row_number + 2)
        def gtf = resolve_manifest_path(row.gtf, manifest_base, 'gene_model' in selected, row_number + 2)
        def sj = resolve_manifest_path(row.sj_out_tab, manifest_base, 'splice_junction' in selected && !derive_junctions_enabled(), row_number + 2)
        def supplied_bai = resolve_manifest_path(row.bai, manifest_base, false, row_number + 2)
        def bai = resolve_current_bam_index(bam, supplied_bai, row.entity_id)
        def meta = [id: row.entity_id, safe_id: safe, entity_type: row.entity_type,
                    gca_accession: row.gca_accession, assembly_release: row.assembly_release,
                    source_entity_ids: row.source_entity_ids ?: '', run_accession: row.run_accession ?: '',
                    sample_id: row.sample_id ?: '', track_types: selected,
                    derive_junctions: derive_junctions_enabled(), contains_lts_path: contains_lts_path]
        records << tuple(meta, bam, gtf, sj, bai)
    }
    records
}

workflow {
    validateParameters()
    if (!params.input_manifest) error '--input_manifest is required'
    if (!params.chrom_sizes) error '--chrom_sizes is required'
    if (!params.assembly_release) error '--assembly_release is required'
    def base_dir = params.manifest_base_dir ? file(params.manifest_base_dir) : new File(params.input_manifest.toString()).getAbsoluteFile().getParentFile()
    if (!base_dir.exists()) error "Manifest base directory does not exist: ${base_dir}"

    def entity_records = manifest_entities(params.input_manifest, base_dir, params.track_types == 'all' ? null : params.track_types)
    def has_lts_inputs = entity_records.any { meta, _bam, _gtf, _sj, _bai -> meta.contains_lts_path }
    if (has_lts_inputs) {
        if (!params.lts_stage_dir) error '--lts_stage_dir is required when any manifest input is under /lts'
        def stage_dir_path = new File(params.lts_stage_dir.toString()).toPath().toAbsolutePath().normalize().toString()
        if (!(stage_dir_path.startsWith('/hps/') || stage_dir_path.startsWith('/nfs/')))
            error "--lts_stage_dir must be under /hps or /nfs when staging /lts inputs: ${stage_dir_path}"
    }

    source_entities = channel.fromList(entity_records)
    STAGE_LTS_INPUTS(source_entities)
    entities = STAGE_LTS_INPUTS.out.entities.map { meta, bam, gtf, sj, bai ->
        tuple(meta, bam, gtf, sj, resolve_current_bam_index(bam, bai, meta.id, false))
    }
    chrom_sizes = file(params.chrom_sizes, checkIfExists: true)
    PREFLIGHT_TRACK_INPUTS(chrom_sizes, params.assembly_release)

    bam_index_entities = entities
        .filter { meta, bam, gtf, sj, bai ->
            bam && !bai && ('coverage' in meta.track_types || (meta.derive_junctions && 'splice_junction' in meta.track_types))
        }
        .map { meta, bam, gtf, sj, bai -> tuple(meta, bam) }
    generated_indexes = INDEX_BAM(bam_index_entities)

    supplied_index_entities = entities
        .filter { meta, bam, gtf, sj, bai ->
            bam && bai && ('coverage' in meta.track_types || (meta.derive_junctions && 'splice_junction' in meta.track_types))
        }
        .map { meta, bam, gtf, sj, bai -> tuple(meta, bam, bai) }
    indexed_entities = supplied_index_entities.mix(generated_indexes.indexed)

    coverage_entities = indexed_entities
        .filter { meta, bam, bai -> 'coverage' in meta.track_types }
        .map { meta, bam, bai -> tuple(meta, bam, bai, chrom_sizes) }
    gene_entities = entities.filter { meta, bam, gtf, sj, bai -> 'gene_model' in meta.track_types }
        .map { meta, bam, gtf, sj, bai -> tuple(meta, gtf, chrom_sizes) }
    splice_entities = entities.filter { meta, bam, gtf, sj, bai -> 'splice_junction' in meta.track_types }
        .map { meta, bam, gtf, sj, bai -> tuple(meta, derive_junctions_enabled() ? bam : sj, bam ?: sj, chrom_sizes) }

    coverage = PREPARE_COVERAGE_TRACKS(coverage_entities)
    genes = PREPARE_GENE_MODEL_TRACKS(gene_entities)
    splice = PREPARE_SPLICE_TRACKS(splice_entities)

    result_channels = coverage.results.map { meta, artifact, result, provenance, exception -> result }
        .mix(genes.results.map { meta, artifact, result, provenance, exception -> result })
        .mix(splice.results.map { meta, artifact, result, provenance, exception -> result })
    version_channels = PREFLIGHT_TRACK_INPUTS.out.versions
        .mix(STAGE_LTS_INPUTS.out.versions)
        .mix(coverage.versions)
        .mix(generated_indexes.versions)
        .mix(genes.versions)
        .mix(splice.versions)
    ASSEMBLE_TRACK_SETS(result_channels, version_channels,
        file("${workflow.launchDir}/pipelines/track_file_preparation/bin/assemble_track_reports.py", checkIfExists: true))
}
