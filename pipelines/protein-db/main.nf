#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'
include { BUILD_PROTEIN_DATABASE } from './subworkflows/build_protein_database.nf'

def parse_boolean(value, field, row_number) {
    def normalised = value?.toString()?.trim()?.toLowerCase()
    if (!normalised || normalised == '-') {
        return false
    }
    if (normalised in ['true', 'yes', '1']) {
        return true
    }
    if (normalised in ['false', 'no', '0']) {
        return false
    }
    error "Invalid '${field}' value '${value}' on build manifest row ${row_number}; use true or false"
}

def manifest_value(values, column_index, column) {
    def position = column_index[column]
    position == null || position >= values.size() ? '' : values[position].trim()
}

def load_build_manifest(manifest_file) {
    def lines = manifest_file.readLines().findAll { line ->
        def stripped = line.trim()
        stripped && !stripped.startsWith('#')
    }
    if (lines.size() < 2) {
        error "Build manifest ${manifest_file} must contain a header and at least one build row"
    }

    def header = lines[0].split('\\t', -1).collect { column -> column.trim() }
    def required = ['build_id', 'source', 'scope', 'taxon_id']
    def missing = required.findAll { column -> !header.contains(column) }
    if (missing) {
        error "Build manifest is missing required columns: ${missing.join(', ')}"
    }

    def index = [:]
    header.eachWithIndex { column, position -> index[column] = position }
    def build_ids = [] as Set
    def records = []

    lines.drop(1).eachWithIndex { line, offset ->
        def row_number = offset + 2
        def values = line.split('\\t', -1)
        def build_id = manifest_value(values, index, 'build_id')
        def source = manifest_value(values, index, 'source').toLowerCase()
        def scope = manifest_value(values, index, 'scope').toLowerCase()
        def taxon_id = manifest_value(values, index, 'taxon_id')
        def taxon_name = manifest_value(values, index, 'taxon_name') ?: build_id
        def query_extra = manifest_value(values, index, 'query_extra')
        def include_isoforms = parse_boolean(manifest_value(values, index, 'include_isoforms'), 'include_isoforms', row_number)

        if (!build_id || !(build_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "Invalid or missing build_id on build manifest row ${row_number}: '${build_id}'"
        }
        if (!build_ids.add(build_id)) {
            error "Duplicate build_id '${build_id}' on build manifest row ${row_number}"
        }
        if (!(source in ['uniprot', 'swissprot'])) {
            error "Invalid source '${source}' on build manifest row ${row_number}; use uniprot or swissprot"
        }
        if (!(scope in ['species', 'clade'])) {
            error "Invalid scope '${scope}' on build manifest row ${row_number}; use species or clade"
        }
        if (!(taxon_id ==~ /[0-9]+/)) {
            error "Invalid taxon_id '${taxon_id}' on build manifest row ${row_number}; use an NCBI taxonomy ID"
        }
        if (query_extra.contains('\\t') || query_extra.contains('\\n') || query_extra.contains('\\r')) {
            error "query_extra on build manifest row ${row_number} must not contain tabs or newlines"
        }

        def meta = [
            id: build_id,
            build_id: build_id,
            source: source,
            scope: scope,
            taxon_id: taxon_id,
            taxon_name: taxon_name,
            query_extra: query_extra == '-' ? '' : query_extra,
            include_isoforms: include_isoforms
        ]
        records << tuple(meta, meta)
    }
    records
}

workflow {
    validateParameters()

    def manifest_file = file(params.build_manifest, checkIfExists: true)
    def fetch_script = file("${projectDir}/bin/fetch_uniprot_fasta.py", checkIfExists: true)
    def builds = load_build_manifest(manifest_file)
    build_ch = channel.fromList(builds)

    BUILD_PROTEIN_DATABASE(
        build_ch,
        params.uniprot_api,
        params.page_size,
        fetch_script
    )
}
