#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { DOWNLOAD_FASTQ } from './modules/download_fastq'

def safe_id(value, field, row_number) {
    def result = (value ?: '').toString().trim()
    if (!(result ==~ /[A-Za-z0-9_.-]+/))
        error "manifest row ${row_number}: ${field} must contain only letters, numbers, dot, underscore, or hyphen"
    result
}

def manifest_path(value, root, field, row_number) {
    def raw = (value ?: '').toString().trim()
    if (!raw) error "manifest row ${row_number}: ${field} is required"
    // Accept local files, explicit URIs, and the bare ENA FTP host paths used
    // by the project batch tables. Remote URIs are acquired by fastq-dl later.
    if (raw.startsWith('ftp.')) raw = "https://${raw}"
    def remote = raw ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/
    if (remote) return raw
    def candidate = raw.startsWith('/') ? file(raw) : file("${root}/${raw}")
    if (!candidate.exists()) error "manifest row ${row_number}: ${field} does not exist: ${candidate}"
    candidate.toAbsolutePath().normalize()
}

def reference_cache_root() {
    params.reference_cache_dir ?: "${params.outdir}/per_assembly"
}

def reference_cache_assets(meta) {
    def directory = file("${reference_cache_root()}/${meta.species_id}/${meta.assembly_accession}/reference")
    def fasta = file("${directory}/reference.fa")
    if (!fasta.exists()) {
        def legacy_fastas = directory.toFile().listFiles()?.findAll { asset -> asset.isFile() && asset.name.endsWith('.fa') }
        if (legacy_fastas?.size() == 1) fasta = legacy_fastas[0].toPath()
    }
    def minimap_index = file("${directory}/reference.mmi")
    def fasta_index = file("${directory}/reference.fa.fai")
    def version = file("${directory}/reference_cache.version")
    def assets_present = [fasta, minimap_index, fasta_index].every { asset ->
        asset.exists() && asset.toFile().length() > 0
    }
    def versioned_cache = version.exists() && version.text.trim() == 'reference-cache-v1'
    def legacy_cache = !version.exists() && [
        file("${directory}/fetch_genome.versions.yml"),
        file("${directory}/minimap2_index.versions.yml"),
        file("${directory}/samtools_faidx.versions.yml")
    ].every { asset -> asset.exists() && asset.toFile().length() > 0 }
    def complete = assets_present && (versioned_cache || legacy_cache)
    [complete: complete, fasta: fasta, minimap_index: minimap_index, fasta_index: fasta_index]
}

process FETCH_GENOME {
    label 'genome_download'
    tag "${meta.species_id}/${meta.assembly_id}"
    container 'docker.io/library/python@sha256:17b9724493d916f56565c131f389faff2597f46dec6ce09516f95125bc9f519f'
    publishDir path: { "${reference_cache_root()}/${meta.species_id}/${meta.assembly_accession}/reference" }, mode: 'copy', saveAs: { filename ->
        filename.endsWith('.fa') ? 'reference.fa' : filename
    }

    input:
    tuple val(key), val(meta), val(assembly_accession), path(fetch_script, stageAs: 'fetch_genome.py')

    output:
    tuple val(key), val(meta), path('*.fa'), emit: reference
    path '*_assembly_report.txt', optional: true, emit: assembly_report
    path 'reference_cache.version', emit: cache_version
    path 'fetch_genome.versions.yml', emit: versions

    script:
    """
    python3 fetch_genome.py --gca '${assembly_accession}' --output_dir .
    printf 'reference-cache-v1\\n' > reference_cache.version
    printf 'fetch_genome.py: %s\\n' "\$(python3 fetch_genome.py --version)" > fetch_genome.versions.yml
    """

    stub:
    """
    printf '>chr1\\nACGT\\n' > '${assembly_accession.replace('.', 'v')}.fa'
    printf 'reference-cache-v1\\n' > reference_cache.version
    printf 'stub\\n' > fetch_genome.versions.yml
    """
}

process BUILD_MINIMAP2_INDEX {
    label 'index'
    tag "${meta.species_id}/${meta.assembly_id}"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:9352893c7f30ec6a3d6e5591d14a4034260a9489537d6577a0a23ae8ee8fd259'
    publishDir path: { "${reference_cache_root()}/${meta.species_id}/${meta.assembly_accession}/reference" }, mode: 'copy', saveAs: { filename ->
        filename in ['reference.mmi', 'minimap2_index.versions.yml'] ? filename : null
    }

    input:
    tuple val(key), val(meta), path(reference, stageAs: 'reference.fa')

    output:
    tuple val(key), val(meta), path('reference.mmi'), emit: index
    path 'minimap2_index.versions.yml', emit: versions

    script:
    """
    minimap2 -x splice -d reference.mmi reference.fa
    cat > minimap2_index.versions.yml <<VERSIONS
"${task.process}":
    minimap2: \$(minimap2 --version)
VERSIONS
    """

    stub:
    """
    printf 'stub index\\n' > reference.mmi
    touch minimap2_index.versions.yml
    """
}

process INDEX_REFERENCE_FASTA {
    label 'fasta_index'
    tag "${meta.species_id}/${meta.assembly_id}"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:9352893c7f30ec6a3d6e5591d14a4034260a9489537d6577a0a23ae8ee8fd259'
    publishDir path: { "${reference_cache_root()}/${meta.species_id}/${meta.assembly_accession}/reference" }, mode: 'copy', saveAs: { filename ->
        filename in ['reference.fa.fai', 'samtools_faidx.versions.yml'] ? filename : null
    }

    input:
    tuple val(key), val(meta), path(reference, stageAs: 'reference.fa')

    output:
    tuple val(key), val(meta), path('reference.fa.fai'), emit: index
    path 'samtools_faidx.versions.yml', emit: versions

    script:
    """
    samtools faidx reference.fa
    cat > samtools_faidx.versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    printf 'chr1\\t1000\\t0\\t1000\\t1001\\n' > reference.fa.fai
    touch samtools_faidx.versions.yml
    """
}

process ARCHIVE_INPUT_MANIFEST {
    label 'manifest'
    tag 'input_manifest'
    publishDir path: { "${params.outdir}/pipeline_info" }, mode: 'copy'

    input:
    path manifest

    output:
    path 'input_manifest.tsv'
    path 'input_manifest.tsv.sha256'

    script:
    """
    cp ${manifest} input_manifest.tsv
    sha256sum input_manifest.tsv > input_manifest.tsv.sha256
    """

    stub:
    """
    cp ${manifest} input_manifest.tsv
    printf 'stub  input_manifest.tsv\\n' > input_manifest.tsv.sha256
    """
}

process ALIGN_AND_SORT_BAM {
    label 'alignment'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:9352893c7f30ec6a3d6e5591d14a4034260a9489537d6577a0a23ae8ee8fd259'
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/alignment" }, mode: 'copy', saveAs: { filename ->
        filename in ['aligned.sorted.bam', 'alignment.versions.yml'] ? filename : null
    }

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('aligned.sorted.bam'), emit: aligned
    path 'alignment.versions.yml', emit: versions

    script:
    // minimap2 uses an extra I/O thread; reserve the remaining two slots for samtools sort.
    def threads = Math.max(1, task.cpus - 3)
    """
    set -o pipefail
    minimap2 -t ${threads} -ax ${meta.minimap2_preset} --secondary=no ${index} ${reads} \\
      | samtools sort -@ 1 -o aligned.sorted.bam -
    cat > alignment.versions.yml <<VERSIONS
"${task.process}":
    minimap2: \$(minimap2 --version)
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    printf 'BAM\\n' > aligned.sorted.bam
    touch alignment.versions.yml
    """
}

process INDEX_ALIGNMENT_BAM {
    label 'bam_index'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:9352893c7f30ec6a3d6e5591d14a4034260a9489537d6577a0a23ae8ee8fd259'
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/alignment" }, mode: 'copy', saveAs: { filename ->
        filename in ['aligned.sorted.bam.csi', 'bam_index.versions.yml'] ? filename : null
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('aligned.sorted.bam.csi'), emit: index
    path 'bam_index.versions.yml', emit: versions

    script:
    """
    samtools index -c -o aligned.sorted.bam.csi ${bam}
    cat > bam_index.versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    printf 'CSI\\n' > aligned.sorted.bam.csi
    touch bam_index.versions.yml
    """
}

process QC_ALIGNMENT_BAM {
    label 'alignment_qc'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:9352893c7f30ec6a3d6e5591d14a4034260a9489537d6577a0a23ae8ee8fd259'
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/alignment" }, mode: 'copy', saveAs: { filename ->
        filename in ['flagstat.txt', 'idxstats.tsv', 'alignment_qc.versions.yml'] ? filename : null
    }

    input:
    tuple val(meta), path(bam), path(reference_fai, stageAs: 'reference.fa.fai')

    output:
    path 'flagstat.txt', emit: flagstat
    path 'idxstats.tsv', emit: idxstats
    path 'alignment_qc.versions.yml', emit: versions

    script:
    """
    samtools quickcheck -v ${bam}
    samtools idxstats ${bam} > idxstats.tsv
    cut -f1,2 reference.fa.fai | sort > expected_contigs.tsv
    awk '\$1 != "*" {print \$1 "\\t" \$2}' idxstats.tsv | sort > observed_contigs.tsv
    diff -u expected_contigs.tsv observed_contigs.tsv
    samtools flagstat ${bam} > flagstat.txt
    cat > alignment_qc.versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    printf 'chr1\\t10\\t0\\t0\\t0\\n*\\t0\\t0\\t0\\t0\\n' > idxstats.tsv
    printf 'expected\\n' > flagstat.txt
    touch alignment_qc.versions.yml
    """
}

process STRINGTIE3_ASSEMBLE {
    label 'stringtie_assembly'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34d25ca28cb224a73658053a1f31'
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/stringtie3/assembly" }, mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.run_id}.assembled.gtf"), path('assembly_status.txt'), emit: assembled
    path 'gene_abundances.tsv', emit: gene_abundances
    path 'versions.yml', emit: versions

    script:
    def nascent_args = meta.nascent_mode == 'mature' ? '-N' : (meta.nascent_mode == 'retain' ? '--nasc' : '')
    """
    stringtie -L ${nascent_args} -p ${task.cpus} \\
      -A gene_abundances.tsv \\
      -o '${meta.run_id}.assembled.gtf' ${bam}
    if awk -F '\t' '\$3 == "transcript" { found=1; exit } END { exit !found }' '${meta.run_id}.assembled.gtf'; then
      printf 'TRANSCRIPTS_FOUND\\n' > assembly_status.txt
    else
      printf 'NO_TRANSCRIPTS\\n' > assembly_status.txt
    fi
    printf 'stringtie: %s\\n' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """

    stub:
    """
    printf 'chr1\\tStringTie\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' > '${meta.run_id}.assembled.gtf'
    printf 'TRANSCRIPTS_FOUND\\n' > assembly_status.txt
    printf 'gene\\tcoverage\\tFPKM\\tTPM\\ng1\\t1\\t1\\t1\\n' > gene_abundances.tsv
    touch versions.yml
    """
}

process STRINGTIE3_MERGE {
    label 'stringtie_merge'
    tag "${meta.species_id}/${meta.assembly_id}"
    container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34d25ca28cb224a73658053a1f31'
    publishDir path: { "${params.outdir}/per_assembly/${meta.species_id}/${meta.assembly_id}/stringtie3/merge" }, mode: 'copy'

    input:
    tuple val(meta), path(gtfs)

    output:
    tuple val(meta), path('merged.gtf'), emit: merged
    path 'versions.yml', emit: versions

    script:
    """
    stringtie --merge -p ${task.cpus} -o merged.gtf ${gtfs}
    printf 'stringtie: %s\\n' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """

    stub:
    """
    printf 'chr1\\tStringTie\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' > merged.gtf
    touch versions.yml
    """
}

process STRINGTIE3_QUANTIFY {
    label 'stringtie_quantification'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34d25ca28cb224a73658053a1f31'
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/stringtie3/quantification" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(merged_gtf)

    output:
    tuple val(meta), path('quantified.gtf'), emit: quantified
    path 'gene_abundances.tsv', emit: gene_abundances
    path 'ballgown/*.ctab', emit: ballgown_tables
    path 'versions.yml', emit: versions

    script:
    """
    stringtie -L -e -p ${task.cpus} -G ${merged_gtf} \\
      -A gene_abundances.tsv -b ballgown -o quantified.gtf ${bam}
    printf 'stringtie: %s\\n' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """

    stub:
    """
    mkdir -p ballgown
    printf 'chr1\\tStringTie\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' > quantified.gtf
    printf 'gene\\tcoverage\\tFPKM\\tTPM\\ng1\\t1\\t1\\t1\\n' > gene_abundances.tsv
    for table in e2t.ctab e_data.ctab i2t.ctab i_data.ctab t_data.ctab; do printf 'stub\\n' > "ballgown/\$table"; done
    touch versions.yml
    """
}

workflow {
    validateParameters()
    def defer_merge = params.defer_merge.toString().trim().toLowerCase() == 'true'
    if (!(params.run_mode in ['batch', 'merge']))
        error "run_mode must be 'batch' or 'merge'"
    if (params.run_mode == 'merge' && defer_merge)
        error 'defer_merge is only valid with run_mode=batch'
    log.info(paramsSummaryLog(workflow))

    manifest_file = file(params.manifest, checkIfExists: true).toAbsolutePath().normalize()
    ARCHIVE_INPUT_MANIFEST(channel.fromPath(manifest_file))
    merge_gtfs = Channel.empty()
    quantify_bams = Channel.empty()
    if (params.run_mode == 'merge') {
        merge_manifest_rows = channel.fromPath(manifest_file)
            .splitCsv(header: true, sep: '\t')
            .collect()
            .map { rows ->
                if (!rows) error 'merge manifest has no data rows'
                def run_ids = [] as Set
                def assembly_accessions = [:]
                rows.eachWithIndex { row, idx ->
                    def row_number = idx + 2
                    def species_id = safe_id(row.species_id, 'species_id', row_number)
                    def assembly_id = safe_id(row.assembly_id, 'assembly_id', row_number)
                    def run_id = safe_id(row.run_id, 'run_id', row_number)
                    if (!run_ids.add(run_id)) error "merge manifest contains duplicate run_id '${run_id}'"
                    def assembly_accession = (row.assembly_accession ?: '').toString().trim()
                    if (!(assembly_accession ==~ /GCA_\d{9}\.\d+/))
                        error "merge manifest row ${row_number}: assembly_accession must be a GCA accession"
                    def assembly_key = "${species_id}::${assembly_id}"
                    if (assembly_accessions.containsKey(assembly_key) && assembly_accessions[assembly_key] != assembly_accession)
                        error "merge manifest rows for ${assembly_key} use different assembly_accessions"
                    assembly_accessions[assembly_key] = assembly_accession
                    def assembly_status = (row.assembly_status ?: '').toString().trim()
                    if (!(assembly_status in ['TRANSCRIPTS_FOUND', 'NO_TRANSCRIPTS', 'ASSEMBLY_MISSING']))
                        error "merge manifest row ${row_number}: assembly_status must be TRANSCRIPTS_FOUND, NO_TRANSCRIPTS, or ASSEMBLY_MISSING"
                    def bam_path = (row.bam_path ?: '').toString().trim()
                    if (!bam_path) error "merge manifest row ${row_number}: bam_path is required"
                    row._bam = file(bam_path, checkIfExists: true)
                    def gtf_path = (row.gtf_path ?: '').toString().trim()
                    if (assembly_status != 'ASSEMBLY_MISSING' && !gtf_path)
                        error "merge manifest row ${row_number}: gtf_path is required for ${assembly_status}"
                    row._gtf = gtf_path ? file(gtf_path, checkIfExists: true) : null
                    row._normalized = [
                        key: assembly_key,
                        species_id: species_id,
                        assembly_id: assembly_id,
                        assembly_accession: assembly_accession,
                        run_id: run_id,
                        sample_id: safe_id(row.sample_id ?: run_id, 'sample_id', row_number)
                    ]
                    row._assembly_status = assembly_status
                }
                rows
            }

        merge_gtfs = merge_manifest_rows.flatMap { rows ->
            rows.findAll { it._assembly_status == 'TRANSCRIPTS_FOUND' }
                .collect { row -> tuple(row._normalized.key, row._normalized, row._gtf) }
        }
        quantify_bams = merge_manifest_rows.flatMap { rows ->
            rows.collect { row -> tuple(row._normalized.key, row._normalized, row._bam) }
        }
    } else {
      manifest_rows = channel.fromPath(manifest_file)
        .splitCsv(header: true, sep: '\t')
        .collect()
        .map { rows ->
            if (!rows) error 'manifest has no data rows'
            def root = manifest_file.parent
            def run_ids = [] as Set
            def assembly_accessions = [:]
            def assembly_nascent_modes = [:]
            rows.eachWithIndex { row, idx ->
                def row_number = idx + 2
                def species_id = safe_id(row.species_id, 'species_id', row_number)
                def assembly_id = safe_id(row.assembly_id, 'assembly_id', row_number)
                def run_id = safe_id(row.run_id, 'run_id', row_number)
                if (!run_ids.add(run_id)) error "manifest contains duplicate run_id '${run_id}'"
                def technology = (row.technology ?: '').toString().trim().toLowerCase()
                if (!(technology in ['ont', 'pacbio'])) error "manifest row ${row_number}: technology must be ont or pacbio"
                def preset = row.minimap2_preset?.toString()?.trim() ?: (technology == 'ont' ? 'splice' : 'splice:hq')
                if (!(preset in ['splice', 'splice:hq'])) error "manifest row ${row_number}: minimap2_preset must be splice or splice:hq"
                def nascent_mode = row.nascent_mode?.toString()?.trim()?.toLowerCase() ?: 'off'
                if (!(nascent_mode in ['off', 'mature', 'retain'])) error "manifest row ${row_number}: nascent_mode must be off, mature, or retain"
                def reads = manifest_path(row.fastq, root, 'fastq', row_number)
                if (reads instanceof CharSequence && !(run_id ==~ /(?:ERR|SRR|DRR)\d+/))
                    error "manifest row ${row_number}: remote FASTQ acquisition requires an ENA/SRA run accession in run_id"
                def assembly_accession = (row.assembly_accession ?: '').toString().trim()
                if (!(assembly_accession ==~ /GCA_\d{9}\.\d+/))
                    error "manifest row ${row_number}: assembly_accession must be a GCA accession such as GCA_000001405.29"
                def assembly_key = "${species_id}::${assembly_id}"
                if (assembly_accessions.containsKey(assembly_key) && assembly_accessions[assembly_key] != assembly_accession)
                    error "manifest rows for ${assembly_key} must use the same assembly_accession"
                if (assembly_nascent_modes.containsKey(assembly_key) && assembly_nascent_modes[assembly_key] != nascent_mode)
                    error "manifest rows for ${assembly_key} must use the same nascent_mode so their transcript sets can be merged consistently"
                assembly_accessions[assembly_key] = assembly_accession
                assembly_nascent_modes[assembly_key] = nascent_mode
                row._reads = reads
                row._assembly_accession = assembly_accession
                row._normalized = [
                    key: assembly_key,
                    species_id: species_id,
                    assembly_id: assembly_id,
                    assembly_accession: assembly_accession,
                    run_id: run_id,
                    sample_id: safe_id(row.sample_id ?: run_id, 'sample_id', row_number),
                    technology: technology,
                    minimap2_preset: preset,
                    nascent_mode: nascent_mode
                ]
            }
            rows.collect { row ->
                def meta = row._normalized
                tuple(meta.key, meta, row._reads, row._assembly_accession)
            }
        }

    reference_specs = manifest_rows.flatMap { records ->
        records.groupBy { record -> record[0] }.collect { key, records_for_assembly ->
            def first = records_for_assembly[0]
            def first_meta = first[1]
            def reference_meta = [
                key: key,
                species_id: first_meta.species_id,
                assembly_id: first_meta.assembly_id,
                assembly_accession: first[3]
            ]
            tuple(key, reference_meta, first[3], file("${projectDir}/bin/fetch_genome.py"))
        }
    }
    reference_inputs = reference_specs
        .filter { key, meta, accession, fetch_script -> !reference_cache_assets(meta).complete }
    cached_references = reference_specs
        .filter { key, meta, accession, fetch_script -> reference_cache_assets(meta).complete }
        .map { key, meta, accession, fetch_script ->
            def assets = reference_cache_assets(meta)
            tuple(key, meta, assets.fasta, assets.minimap_index, assets.fasta_index)
        }

    FETCH_GENOME(reference_inputs)
    fetched_reference_sequences = FETCH_GENOME.out.reference
    BUILD_MINIMAP2_INDEX(fetched_reference_sequences.map { key, meta, reference -> tuple(key, meta, reference) })
    INDEX_REFERENCE_FASTA(fetched_reference_sequences.map { key, meta, reference -> tuple(key, meta, reference) })

    minimap_indexes = BUILD_MINIMAP2_INDEX.out.index.mix(
        cached_references.map { key, meta, _reference, minimap_index, _fasta_index -> tuple(key, meta, minimap_index) }
    )
    fasta_indexes = INDEX_REFERENCE_FASTA.out.index.mix(
        cached_references.map { key, meta, _reference, _minimap_index, fasta_index -> tuple(key, meta, fasta_index) }
    )

    run_inputs = manifest_rows
        .flatMap { records -> records }
    local_read_inputs = run_inputs
        .filter { key, meta, reads, _assembly_accession -> reads instanceof java.nio.file.Path }
        .map { key, meta, reads, _assembly_accession -> tuple(meta, reads) }
    remote_read_inputs = run_inputs
        .filter { key, meta, reads, _assembly_accession -> reads instanceof CharSequence }
        .map { key, meta, reads, _assembly_accession -> tuple(meta, reads.toString()) }
    DOWNLOAD_FASTQ(remote_read_inputs)
    read_inputs = local_read_inputs.mix(DOWNLOAD_FASTQ.out.reads)

    alignment_inputs = read_inputs
        .map { meta, reads -> tuple(meta.key, meta, reads) }
        .combine(minimap_indexes, by: 0)
        .map { key, meta, reads, _assembly_meta, index -> tuple(meta, reads, index) }
    ALIGN_AND_SORT_BAM(alignment_inputs)
    INDEX_ALIGNMENT_BAM(ALIGN_AND_SORT_BAM.out.aligned)

    alignment_qc_inputs = ALIGN_AND_SORT_BAM.out.aligned
        .map { meta, bam -> tuple(meta.key, meta, bam) }
        .combine(fasta_indexes.map { key, _meta, reference_fai -> tuple(key, reference_fai) }, by: 0)
        .map { key, meta, bam, reference_fai -> tuple(meta, bam, reference_fai) }
    QC_ALIGNMENT_BAM(alignment_qc_inputs)

    STRINGTIE3_ASSEMBLE(ALIGN_AND_SORT_BAM.out.aligned.map { meta, bam -> tuple(meta, bam) })

    if (!defer_merge) {
        merge_gtfs = STRINGTIE3_ASSEMBLE.out.assembled
            .filter { meta, gtf, status -> status.text.trim() == 'TRANSCRIPTS_FOUND' }
            .map { meta, gtf, status -> tuple(meta.key, meta, gtf) }
        quantify_bams = ALIGN_AND_SORT_BAM.out.aligned
            .map { meta, bam -> tuple(meta.key, meta, bam) }
    } else {
        log.info('Species-level merge and quantification deferred; per-run BAMs, assemblies, and assembly status files are published for a later merge run.')
    }
    }

    if (!defer_merge) {
        merge_groups = merge_gtfs
            .groupTuple(by: 0)
            .map { key, metas, gtfs -> tuple(key, metas[0], gtfs) }
        STRINGTIE3_MERGE(merge_groups.map { key, meta, gtfs -> tuple(meta, gtfs) })

        quantified_inputs = quantify_bams
            .combine(STRINGTIE3_MERGE.out.merged.map { meta, merged_gtf -> tuple(meta.key, merged_gtf) }, by: 0)
            .map { key, meta, bam, merged_gtf -> tuple(meta, bam, merged_gtf) }
        STRINGTIE3_QUANTIFY(quantified_inputs)
    }
}
