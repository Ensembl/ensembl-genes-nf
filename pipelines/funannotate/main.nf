#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'
include { STAGE_FTP_INPUTS } from './subworkflows/stage_ftp_inputs.nf'
include { CONVERT_CRAMS_TO_BAMS } from './modules/convert_crams_to_bams.nf'
include { VALIDATE_GENOME } from './modules/validate_genome.nf'
include { VALIDATE_BAMS } from './modules/validate_bams.nf'
include { INDEX_BAMS } from './modules/index_bams.nf'
include { RUN_FUNANNOTATE } from './modules/run_funannotate.nf'

def split_list(value) {
    if (!value?.toString()?.trim()) return []
    value.toString().split(';', -1).collect { it.trim() }.findAll { it }
}

def safe_id(String value) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def require_ftp(String value, String field, int row) {
    if (!value?.trim()) error "Manifest row ${row} requires ${field}"
    if (!(value ==~ /(?i)ftp:\/\/.+/)) error "Manifest row ${row} ${field} must be an ftp:// URL"
}

def parse_manifest(manifest_file) {
    def lines = manifest_file.readLines()
    if (!lines) error "Input manifest is empty: ${manifest_file}"
    def header = lines[0].split('\t', -1).collect { it.trim() }
    def required = ['cohort_id', 'assembly_accession', 'assembly_release', 'species', 'genome_path', 'cram_ftp_urls']
    def missing = required.findAll { !(it in header) }
    if (missing) error "Manifest is missing required columns: ${missing.join(', ')}"
    if (header.size() != header.toSet().size()) error 'Manifest contains duplicate column names'
    def index = header.withIndex().collectEntries { name, i -> [(name): i] }
    def seen = [] as Set
    def records = []
    lines.drop(1).eachWithIndex { line, offset ->
        def row_number = offset + 2
        if (!line.trim()) return
        def values = line.split('\t', -1)
        if (values.size() != header.size()) error "Manifest row ${row_number} has ${values.size()} fields; expected ${header.size()}"
        def row = header.collectEntries { name -> [(name): values[index[name]].trim()] }
        def key = "${row.assembly_accession}\t${row.cohort_id}"
        if (!seen.add(key)) error "Duplicate cohort_id plus assembly_accession: ${key}"
        ['cohort_id', 'assembly_accession', 'assembly_release', 'species'].each { field ->
            if (!row[field]) error "Manifest row ${row_number} requires ${field}"
        }
        def genome = file(row.genome_path, checkIfExists: true)
        def cram_urls = split_list(row.cram_ftp_urls)
        if (!cram_urls) error "Manifest row ${row_number} requires at least one CRAM URL"
        cram_urls.each { require_ftp(it, 'cram_ftp_urls', row_number) }
        def crai_urls = split_list(row.crai_ftp_urls)
        if (crai_urls && crai_urls.size() != cram_urls.size()) error "Manifest row ${row_number} crai_ftp_urls must align with cram_ftp_urls"
        def cram_checksums = split_list(row.cram_sha256s)
        if (cram_checksums && cram_checksums.size() != cram_urls.size()) error "Manifest row ${row_number} cram_sha256s must align with cram_ftp_urls"
        cram_checksums.each { if (!(it ==~ /(?i)[0-9a-f]{64}/)) error "Manifest row ${row_number} contains an invalid CRAM SHA-256" }
        def optional_urls = ['transcript_evidence_ftp_urls', 'protein_evidence_ftp_urls']
        optional_urls.each { field -> split_list(row[field]).each { require_ftp(it, field, row_number) } }
        ['pasa_gff_ftp_url', 'parameters_json_ftp_url', 'stringtie_ftp_url', 'protein_alignments_ftp_url', 'trnascan_ftp_url'].each { field -> if (row[field]) require_ftp(row[field], field, row_number) }
        def safe = safe_id(row.cohort_id)
        def files = [
            genome: [path: genome, name: genome.name],
            crams: cram_urls.withIndex().collect { url, i -> [url: url, checksum: cram_checksums ? cram_checksums[i] : '', name: "cram_${String.format('%02d', i + 1)}.cram"] },
            crais: crai_urls.withIndex().collect { url, i -> [url: url, checksum: '', name: "cram_${String.format('%02d', i + 1)}.cram.crai"] },
            transcripts: split_list(row.transcript_evidence_ftp_urls).withIndex().collect { url, i -> [url: url, checksum: '', name: "transcript_${String.format('%02d', i + 1)}.fa"] },
            proteins: split_list(row.protein_evidence_ftp_urls).withIndex().collect { url, i -> [url: url, checksum: '', name: "protein_${String.format('%02d', i + 1)}.fa"] },
            pasa_gff: row.pasa_gff_ftp_url ? [[url: row.pasa_gff_ftp_url, checksum: '', name: 'pasa.gff3']] : [],
            parameters_json: row.parameters_json_ftp_url ? [[url: row.parameters_json_ftp_url, checksum: '', name: 'parameters.json']] : [],
            stringtie: row.stringtie_ftp_url ? [[url: row.stringtie_ftp_url, checksum: '', name: 'stringtie.gtf']] : [],
            protein_alignments: row.protein_alignments_ftp_url ? [[url: row.protein_alignments_ftp_url, checksum: '', name: 'protein_alignments.gff3']] : [],
            trnascan: row.trnascan_ftp_url ? [[url: row.trnascan_ftp_url, checksum: '', name: 'trnascan.gff']] : []
        ]
        records << [
            id: row.cohort_id, safe_id: safe, assembly_accession: row.assembly_accession,
            assembly_release: row.assembly_release, species: row.species, isolate: row.isolate ?: '',
            files: files, genome_name: 'genome' + (genome.name.toLowerCase().endsWith('.gz') ? '.fa.gz' : '.fa'), bam_count: cram_urls.size(),
            genome_path: genome,
            input_row: row_number
        ]
    }
    if (!records) error 'Input manifest contains no data rows'
    records
}

def validate_extra_args(String value) {
    def owned = ['--input', '--out', '--species', '--isolate', '--rna_bam', '--transcript_evidence', '--protein_evidence', '--pasa_gff', '--parameters', '--cpus']
    def tokens = value?.trim() ? value.trim().split(/\s+/) : []
    def bad = tokens.findAll { token -> owned.any { flag -> token == flag || token.startsWith("${flag}=") } }
    if (bad) error "--extra_predict_args cannot override pipeline-owned arguments: ${bad.unique().join(', ')}"
}

workflow {
    validateParameters()
    validate_extra_args(params.extra_predict_args)
    if (params.parameters && params.parameters.toString().startsWith('ftp://'))
        error '--parameters must be a local parameters JSON path; use parameters_json_ftp_url in the manifest for FTP input'
    def records = parse_manifest(file(params.input_manifest, checkIfExists: true))
    source = channel.fromList(records)
    STAGE_FTP_INPUTS(source.map { record -> tuple(record, record.genome_path) }, params.stage_dir)
    VALIDATE_GENOME(STAGE_FTP_INPUTS.out.staged)
    CONVERT_CRAMS_TO_BAMS(VALIDATE_GENOME.out.validated)
    VALIDATE_BAMS(CONVERT_CRAMS_TO_BAMS.out.converted)
    INDEX_BAMS(VALIDATE_BAMS.out.validated)
    RUN_FUNANNOTATE(
        INDEX_BAMS.out.indexed,
        params.funannotate_db ?: '',
        params.cpus,
        params.augustus_species ?: '',
        params.soft_mask,
        params.min_training_models,
        params.busco_db ?: '',
        params.parameters ?: '',
        params.extra_predict_args ?: ''
    )
}
