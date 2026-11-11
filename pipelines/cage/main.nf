#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'
include { PREPARE_REFERENCES_WORKFLOW } from './subworkflows/prepare_references.nf'
include { PROCESS_READS } from './subworkflows/process_reads.nf'
include { BUILD_CTSS } from './subworkflows/build_ctss.nf'
include { BUILD_SPECIES_ATLAS_WORKFLOW } from './subworkflows/build_species_atlas.nf'
include { BUILD_CROSS_SPECIES_WORKFLOW } from './subworkflows/build_cross_species.nf'

def bool_value(value, field, row_number) {
    def normalised = value?.toString()?.trim()?.toLowerCase()
    if (normalised in ['true', 'yes', '1']) return true
    if (normalised in ['false', 'no', '0', '', '-']) return false
    error "Invalid ${field} '${value}' on manifest row ${row_number}"
}

def required_value(row, field, row_number) {
    def value = row[field]?.toString()?.trim()
    if (!value || value == '-') error "Missing ${field} on manifest row ${row_number}"
    value
}

def manifest_path(value) {
    if (!value) return value
    if (value.startsWith('/')) return value
    def cwd_path = new File(value)
    cwd_path.exists() ? value : "${projectDir}/${value}"
}

def load_protocol_profiles(path) {
    def rows = path.readLines().findAll { line -> line.trim() && !line.trim().startsWith('#') }
    if (rows.size() < 2) error "Protocol profile manifest must contain a header and at least one row"
    def header = rows[0].split('\\t', -1)*.trim()
    def required = ['protocol', 'assay_class', 'trim_args', 'alignment_args', 'ctss_mode']
    def missing = required.findAll { column -> !header.contains(column) }
    if (missing) error "Protocol profile manifest is missing columns: ${missing.join(', ')}"
    def profiles = [:]
    rows.drop(1).eachWithIndex { line, offset ->
        def values = line.split('\\t', -1)
        def row = [:]
        header.eachWithIndex { name, index -> row[name] = index < values.size() ? values[index].trim() : '' }
        def protocol = required_value(row, 'protocol', offset + 2).toLowerCase()
        if (profiles.containsKey(protocol)) error "Duplicate protocol profile '${protocol}'"
        profiles[protocol] = row
    }
    profiles
}

def load_sample_manifest(path, profiles) {
    def rows = path.readLines().findAll { line -> line.trim() && !line.trim().startsWith('#') }
    if (rows.size() < 2) error "Sample manifest must contain a header and at least one row"
    def header = rows[0].split('\t', -1)*.trim()
    def required = ['sample_id', 'species', 'assembly', 'assay_class', 'protocol', 'layout', 'accession', 'fastq_1', 'fastq_2']
    def missing = required.findAll { column -> !header.contains(column) }
    if (missing) error "Sample manifest is missing columns: ${missing.join(', ')}"
    def records = []
    def ids = [] as Set
    rows.drop(1).eachWithIndex { line, offset ->
        def row_number = offset + 2
        def values = line.split('\t', -1)
        def row = [:]
        header.eachWithIndex { name, index -> row[name] = index < values.size() ? values[index].trim() : '' }
        def id = required_value(row, 'sample_id', row_number)
        if (!ids.add(id)) error "Duplicate sample_id '${id}' on manifest row ${row_number}"
        def layout = required_value(row, 'layout', row_number).toUpperCase()
        if (!(layout in ['SE', 'PE'])) error "Unsupported layout '${layout}' on manifest row ${row_number}; use SE or PE"
        def assay = required_value(row, 'assay_class', row_number).toLowerCase()
        if (!(assay in ['cage', 'rampage', 'pro-cap', 'net-cage', 'tss-seq'])) error "Unsupported assay_class '${assay}' on manifest row ${row_number}"
        def protocol = required_value(row, 'protocol', row_number).toLowerCase()
        def protocol_assay = [standard_cage: 'cage', nano_cage: 'cage', helicoscope_cage: 'cage', ss_cage: 'cage', rampage_v1: 'rampage', pro_cap_v1: 'pro-cap', net_cage_v1: 'net-cage', tss_seq_v1: 'tss-seq']
        if (protocol_assay[protocol] != assay || !profiles.containsKey(protocol)) error "Protocol '${protocol}' is not valid for assay_class '${assay}' or has no profile on row ${row_number}"
        def profile = profiles[protocol]
        def fastq1 = row.fastq_1 == '-' ? '' : manifest_path(row.fastq_1)
        def fastq2 = row.fastq_2 == '-' ? '' : manifest_path(row.fastq_2)
        if (!row.accession?.trim() && !fastq1) error "Row ${row_number} needs accession or fastq_1"
        if (layout == 'PE' && !fastq2 && row.accession == '-') error "Paired-end row ${row_number} needs fastq_2 or accession"
        def meta = [
            id: id, sample_id: id, species: required_value(row, 'species', row_number),
            assembly: required_value(row, 'assembly', row_number),
            assay_class: assay, protocol: protocol, trim_args: profile.trim_args ?: '', alignment_args: profile.alignment_args ?: '', ctss_mode: profile.ctss_mode ?: '',
            layout: layout, accession: row.accession == '-' ? '' : row.accession,
            fastq_1: fastq1, fastq_2: fastq2,
            tissue: row.tissue ?: '', condition: row.condition ?: '',
            biological_replicate: row.biological_replicate ?: '',
            checksum_1: row.checksum_1 ?: '', checksum_2: row.checksum_2 ?: ''
        ]
        records << tuple(meta, meta)
    }
    records
}

def load_reference_manifest(path) {
    def rows = path.readLines().findAll { line -> line.trim() && !line.trim().startsWith('#') }
    if (rows.size() < 2) error "Reference manifest must contain a header and at least one row"
    def header = rows[0].split('\t', -1)*.trim()
    def required = ['assembly', 'species', 'fasta', 'gtf', 'chrom_sizes']
    def missing = required.findAll { column -> !header.contains(column) }
    if (missing) error "Reference manifest is missing columns: ${missing.join(', ')}"
    def records = []
    def ids = [] as Set
    rows.drop(1).eachWithIndex { line, offset ->
        def row_number = offset + 2
        def values = line.split('\t', -1)
        def row = [:]
        header.eachWithIndex { name, index -> row[name] = index < values.size() ? values[index].trim() : '' }
        def assembly = required_value(row, 'assembly', row_number)
        if (!ids.add(assembly)) error "Duplicate assembly '${assembly}' on reference manifest row ${row_number}"
        def meta = [id: assembly, assembly: assembly, species: required_value(row, 'species', row_number), taxon_id: row.taxon_id ?: '', source: row.source ?: 'user']
        records << tuple(meta, file(manifest_path(required_value(row, 'fasta', row_number)), checkIfExists: true), file(manifest_path(required_value(row, 'gtf', row_number)), checkIfExists: true), file(manifest_path(required_value(row, 'chrom_sizes', row_number)), checkIfExists: true))
    }
    records
}

workflow {
    validateParameters()
    def protocol_profiles = load_protocol_profiles(file(manifest_path(params.protocol_profiles), checkIfExists: true))
    def samples = load_sample_manifest(file(params.sample_manifest, checkIfExists: true), protocol_profiles)
    def references = load_reference_manifest(file(params.reference_manifest, checkIfExists: true))

    ch_samples = channel.fromList(samples)
    ch_references = channel.fromList(references)

    def reference_assemblies = references.collect { record -> record[0].assembly } as Set
    samples.each { record ->
        if (!reference_assemblies.contains(record[0].assembly)) {
            error "No reference manifest row for sample assembly '${record[0].assembly}'"
        }
    }

    PREPARE_REFERENCES_WORKFLOW(ch_references)
    PROCESS_READS(ch_samples, PREPARE_REFERENCES_WORKFLOW.out.index)
    BUILD_CTSS(PROCESS_READS.out.bam)
    BUILD_SPECIES_ATLAS_WORKFLOW(BUILD_CTSS.out.ctss, PREPARE_REFERENCES_WORKFLOW.out.annotation)

    if (params.run_cross_species) {
        def ch_atlas_dirs = BUILD_SPECIES_ATLAS_WORKFLOW.out.species_atlas.map { _meta, atlas -> atlas }
        BUILD_CROSS_SPECIES_WORKFLOW(ch_atlas_dirs, file(manifest_path(params.orthology_manifest), checkIfExists: true))
    }
}
