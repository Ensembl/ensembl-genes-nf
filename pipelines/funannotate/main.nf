#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'
include { RUN_FUNANNOTATE_PIPELINE } from './workflows/run_funannotate_pipeline.nf'

def split_list(value) {
    if (!value?.toString()?.trim()) {
        return []
    }
    value.toString().split(';', -1).collect { item -> item.trim() }.findAll { item -> item }
}

def safe_id(value: String) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def require_remote_url(value: String, field: String, row: Integer) {
    if (!value?.trim()) {
        error("Manifest row ${row} requires ${field}")
    }
    if (!(value ==~ /(?i)(ftp|https?):\/\/.+/)) {
        error("Manifest row ${row} ${field} must be an ftp:// or https:// URL")
    }
}

def local_spec(value: String, field: String, row: Integer) {
    if (!value?.trim()) {
        error("Manifest row ${row} requires ${field}")
    }
    [path: file(value.trim(), checkIfExists: true), name: value.tokenize('/')[-1]]
}

def local_specs(value: String, field: String, row: Integer) {
    split_list(value).collect { path -> local_spec(path, field, row) }
}

def parse_manifest(manifest_file) {
    def lines = manifest_file.readLines()
    if (!lines) {
        error("Input manifest is empty: ${manifest_file}")
    }
    def header = lines[0].split('\t', -1).collect { item -> item.trim() }
    def required = ['cohort_id', 'assembly_accession', 'species', 'genome_path']
    def missing = required.findAll { field -> !(field in header) }
    if (missing) {
        error("Manifest is missing required columns: ${missing.join(', ')}")
    }
    if (header.size() != header.toSet().size()) {
        error('Manifest contains duplicate column names')
    }
    def index = header.withIndex().collectEntries { name, i -> [(name): i] }
    def seen = [] as Set
    def records = []
    lines
        .drop(1)
        .eachWithIndex { line, offset ->
            def row_number = offset + 2
            if (!line.trim()) {
                return null
            }
            def values = line.split('\t', -1)
            if (values.size() != header.size()) {
                error("Manifest row ${row_number} has ${values.size()} fields; expected ${header.size()}")
            }
            def row = header.collectEntries { name -> [(name): values[index[name]].trim()] }
            def key = "${row.assembly_accession}\t${row.cohort_id}"
            if (!seen.add(key)) {
                error("Duplicate cohort_id plus assembly_accession: ${key}")
            }
            ['cohort_id', 'assembly_accession', 'species'].each { field ->
                if (!row[field]) {
                    error("Manifest row ${row_number} requires ${field}")
                }
            }
            def cram_urls = split_list(row.cram_ftp_urls)
            def cram_dir = row.cram_ftp_dir?.trim() ?: ''
            if (cram_dir && !cram_dir.endsWith('/')) {
                cram_dir += '/'
            }
            if (!cram_urls && !cram_dir) {
                error("Manifest row ${row_number} requires cram_ftp_urls or cram_ftp_dir")
            }
            if (cram_urls && cram_dir) {
                error("Manifest row ${row_number} must provide only one of cram_ftp_urls or cram_ftp_dir")
            }
            cram_urls.each { url -> require_remote_url(url, 'cram_ftp_urls', row_number) }
            if (cram_dir) {
                require_remote_url(cram_dir, 'cram_ftp_dir', row_number)
            }
            def cram_checksums = split_list(row.cram_sha256s)
            if (cram_dir && cram_checksums) {
                error("Manifest row ${row_number} cannot use cram_sha256s with cram_ftp_dir")
            }
            if (cram_checksums && cram_checksums.size() != cram_urls.size()) {
                error("Manifest row ${row_number} cram_sha256s must align with cram_ftp_urls")
            }
            cram_checksums.each { checksum ->
                if (!(checksum ==~ /(?i)[0-9a-f]{64}/)) {
                    error("Manifest row ${row_number} contains an invalid CRAM SHA-256")
                }
            }
            def safe = safe_id(row.cohort_id)
            def files = [genome: local_spec(row.genome_path, 'genome_path', row_number), crams: cram_urls.withIndex().collect { url, i -> [url: url, checksum: cram_checksums ? cram_checksums[i] : '', name: "cram_${String.format('%02d', i + 1)}.cram"] }, cram_dir: cram_dir, transcripts: local_specs(row.transcript_evidence_paths, 'transcript_evidence_paths', row_number), stringtie: row.stringtie_gtf_path ? [local_spec(row.stringtie_gtf_path, 'stringtie_gtf_path', row_number)] : [], scallop: row.scallop_gtf_path ? [local_spec(row.scallop_gtf_path, 'scallop_gtf_path', row_number)] : [], protein_alignments: local_specs(row.protein_alignment_paths, 'protein_alignment_paths', row_number), trnascan: row.trnascan_path ? [local_spec(row.trnascan_path, 'trnascan_path', row_number)] : []]
            def local_files = [files.genome] + files.transcripts + files.stringtie + files.scallop + files.protein_alignments + files.trnascan
            records << [id: row.cohort_id, safe_id: safe, assembly_accession: row.assembly_accession, species: row.species, files: files, local_files: local_files.collect { spec -> spec.path }, genome_name: 'genome.fa', cram_count: cram_urls.size(), input_row: row_number]
        }
    if (!records) {
        error('Input manifest contains no data rows')
    }
    records
}

workflow {
    validateParameters()
    def records = parse_manifest(file(params.input_manifest, checkIfExists: true))
    def cpus = params.cpus.toString().toInteger()
    source = channel.fromList(records)
    RUN_FUNANNOTATE_PIPELINE(
        source.map { record -> tuple(record, record.local_files) },
        params.stage_dir,
        cpus,
    )
}
