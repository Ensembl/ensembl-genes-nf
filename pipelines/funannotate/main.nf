#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'
include { STAGE_FTP_INPUTS } from './subworkflows/stage_ftp_inputs.nf'
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
    def required = ['cohort_id', 'assembly_accession', 'assembly_release', 'species', 'genome_ftp_url', 'bam_ftp_urls']
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
        require_ftp(row.genome_ftp_url, 'genome_ftp_url', row_number)
        def bam_urls = split_list(row.bam_ftp_urls)
        if (!bam_urls) error "Manifest row ${row_number} requires at least one BAM URL"
        bam_urls.each { require_ftp(it, 'bam_ftp_urls', row_number) }
        def bai_urls = split_list(row.bai_ftp_urls)
        if (bai_urls && bai_urls.size() != bam_urls.size()) error "Manifest row ${row_number} bai_ftp_urls must align with bam_ftp_urls"
        def bam_checksums = split_list(row.bam_sha256s)
        if (bam_checksums && bam_checksums.size() != bam_urls.size()) error "Manifest row ${row_number} bam_sha256s must align with bam_ftp_urls"
        bam_checksums.each { if (!(it ==~ /(?i)[0-9a-f]{64}/)) error "Manifest row ${row_number} contains an invalid BAM SHA-256" }
        def optional_urls = ['transcript_evidence_ftp_urls', 'protein_evidence_ftp_urls']
        optional_urls.each { field -> split_list(row[field]).each { require_ftp(it, field, row_number) } }
        ['pasa_gff_ftp_url', 'parameters_json_ftp_url', 'stringtie_ftp_url', 'protein_alignments_ftp_url', 'trnascan_ftp_url'].each { field -> if (row[field]) require_ftp(row[field], field, row_number) }
        def safe = safe_id(row.cohort_id)
        def files = [
            genome: [url: row.genome_ftp_url, checksum: row.genome_sha256, name: row.genome_ftp_url.tokenize('/')[-1]],
            bams: bam_urls.withIndex().collect { url, i -> [url: url, checksum: bam_checksums ? bam_checksums[i] : '', name: "bam_${String.format('%02d', i + 1)}.bam"] },
            bais: bai_urls.withIndex().collect { url, i -> [url: url, checksum: '', name: "bam_${String.format('%02d', i + 1)}.bam.bai"] },
            transcripts: split_list(row.transcript_evidence_ftp_urls).withIndex().collect { url, i -> [url: url, checksum: '', name: "transcript_${String.format('%02d', i + 1)}.fa"] },
            proteins: split_list(row.protein_evidence_ftp_urls).withIndex().collect { url, i -> [url: url, checksum: '', name: "protein_${String.format('%02d', i + 1)}.fa"] },
            pasa_gff: row.pasa_gff_ftp_url ? [[url: row.pasa_gff_ftp_url, checksum: '', name: 'pasa.gff3']] : [],
            parameters_json: row.parameters_json_ftp_url ? [[url: row.parameters_json_ftp_url, checksum: '', name: 'parameters.json']] : [],
            stringtie: row.stringtie_ftp_url ? [[url: row.stringtie_ftp_url, checksum: '', name: 'stringtie.gtf']] : [],
            protein_alignments: row.protein_alignments_ftp_url ? [[url: row.protein_alignments_ftp_url, checksum: '', name: 'protein_alignments.gff3']] : [],
            trnascan: row.trnascan_ftp_url ? [[url: row.trnascan_ftp_url, checksum: '', name: 'trnascan.gff']] : []
        ]
        if (files.genome.checksum && !(files.genome.checksum ==~ /(?i)[0-9a-f]{64}/)) error "Manifest row ${row_number} contains an invalid genome SHA-256"
        records << [
            id: row.cohort_id, safe_id: safe, assembly_accession: row.assembly_accession,
            assembly_release: row.assembly_release, species: row.species, isolate: row.isolate ?: '',
            files: files, genome_name: row.genome_ftp_url.toLowerCase().endsWith('.gz') ? 'genome.fa.gz' : 'genome.fa', bam_count: bam_urls.size(),
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
    STAGE_FTP_INPUTS(source, params.stage_dir)
    VALIDATE_GENOME(STAGE_FTP_INPUTS.out.staged)
    VALIDATE_BAMS(VALIDATE_GENOME.out.validated)
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
