#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

def safe_id(value, field, row_number) {
    def result = (value ?: '').toString().trim()
    if (!(result ==~ /[A-Za-z0-9_.-]+/))
        error "manifest row ${row_number}: ${field} must contain only letters, numbers, dot, underscore, or hyphen"
    result
}

def manifest_path(value, root, field, row_number) {
    def raw = (value ?: '').toString().trim()
    if (!raw) error "manifest row ${row_number}: ${field} is required"
    def candidate = raw.startsWith('/') ? file(raw) : file("${root}/${raw}")
    if (!candidate.exists()) error "manifest row ${row_number}: ${field} does not exist: ${candidate}"
    candidate.toAbsolutePath().normalize()
}

process BUILD_REFERENCE_INDEX {
    label 'index'
    tag "${meta.species_id}/${meta.assembly_id}"
    container params.use_containers ? params.tool_containers.alignment : null
    publishDir path: { "${params.outdir}/per_assembly/${meta.species_id}/${meta.assembly_id}/reference" }, mode: 'copy'

    input:
    tuple val(key), val(meta), path(reference, stageAs: 'reference.fa')

    output:
    tuple val(key), val(meta), path('reference.mmi'), path('reference.fa.fai'), emit: index
    path 'versions.yml', emit: versions

    script:
    """
    minimap2 -x splice -d reference.mmi reference.fa
    samtools faidx reference.fa
    cat > versions.yml <<VERSIONS
"\${task.process}":
    minimap2: \$(minimap2 --version)
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    printf 'stub index\\n' > reference.mmi
    printf 'chr1\\t1000\\t0\\t1000\\t1001\\n' > reference.fa.fai
    touch versions.yml
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

process ALIGN_SORT_INDEX_QC {
    label 'alignment'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container params.use_containers ? params.tool_containers.alignment : null
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/alignment" }, mode: 'copy', saveAs: { filename ->
        filename in ['aligned.sorted.bam', 'aligned.sorted.bam.csi', 'flagstat.txt', 'idxstats.tsv', 'versions.yml'] ? filename : null
    }

    input:
    tuple val(meta), path(reads), path(index), path(reference_fai, stageAs: 'reference.fa.fai')

    output:
    tuple val(meta), path('aligned.sorted.bam'), path('aligned.sorted.bam.csi'), emit: aligned
    path 'versions.yml', emit: versions
    path 'flagstat.txt', emit: flagstat
    path 'idxstats.tsv', emit: idxstats

    script:
    // minimap2 uses an extra I/O thread; reserve the remaining two slots for samtools sort.
    def threads = Math.max(1, task.cpus - 3)
    """
    set -o pipefail
    minimap2 -t ${threads} -ax ${meta.minimap2_preset} --secondary=no ${index} ${reads} \\
      | samtools sort -@ 1 -o aligned.sorted.bam -
    samtools index -c -o aligned.sorted.bam.csi aligned.sorted.bam
    samtools quickcheck -v aligned.sorted.bam
    samtools idxstats aligned.sorted.bam > idxstats.tsv
    cut -f1,2 reference.fa.fai | sort > expected_contigs.tsv
    awk '\$1 != "*" {print \$1 "\\t" \$2}' idxstats.tsv | sort > observed_contigs.tsv
    diff -u expected_contigs.tsv observed_contigs.tsv
    samtools flagstat aligned.sorted.bam > flagstat.txt
    cat > versions.yml <<VERSIONS
"\${task.process}":
    minimap2: \$(minimap2 --version)
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    printf 'BAM\\n' > aligned.sorted.bam
    printf 'CSI\\n' > aligned.sorted.bam.csi
    printf 'chr1\\t10\\t0\\t0\\t0\\n*\\t0\\t0\\t0\\t0\\n' > idxstats.tsv
    printf 'expected\\n' > flagstat.txt
    touch versions.yml
    """
}

process STRINGTIE3_ASSEMBLE {
    label 'stringtie_assembly'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container params.use_containers ? params.tool_containers.stringtie3 : null
    publishDir path: { "${params.outdir}/per_run/${meta.species_id}/${meta.assembly_id}/${meta.run_id}/stringtie3/assembly" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(annotation, stageAs: 'reference.gtf')

    output:
    tuple val(meta), path("${meta.run_id}.assembled.gtf"), emit: assembled
    path 'gene_abundances.tsv', emit: gene_abundances
    path 'fully_covered_reference.gtf', emit: fully_covered
    path 'versions.yml', emit: versions

    script:
    """
    stringtie -L -p ${task.cpus} -G ${annotation} \\
      -A gene_abundances.tsv -C fully_covered_reference.gtf \\
      -o '${meta.run_id}.assembled.gtf' ${bam}
    printf 'stringtie: %s\\n' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """

    stub:
    """
    printf 'chr1\\tStringTie\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' > '${meta.run_id}.assembled.gtf'
    printf 'gene\\tcoverage\\tFPKM\\tTPM\\ng1\\t1\\t1\\t1\\n' > gene_abundances.tsv
    printf 'chr1\\tStringTie\\ttranscript\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' > fully_covered_reference.gtf
    touch versions.yml
    """
}

process STRINGTIE3_MERGE {
    label 'stringtie_merge'
    tag "${meta.species_id}/${meta.assembly_id}"
    container params.use_containers ? params.tool_containers.stringtie3 : null
    publishDir path: { "${params.outdir}/per_assembly/${meta.species_id}/${meta.assembly_id}/stringtie3/merge" }, mode: 'copy'

    input:
    tuple val(meta), path(gtfs), path(annotation, stageAs: 'reference.gtf')

    output:
    tuple val(meta), path('merged.gtf'), emit: merged
    path 'versions.yml', emit: versions

    script:
    """
    stringtie --merge -p ${task.cpus} -G ${annotation} -o merged.gtf ${gtfs}
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
    container params.use_containers ? params.tool_containers.stringtie3 : null
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
    log.info(paramsSummaryLog(workflow))

    manifest_file = file(params.manifest, checkIfExists: true).toAbsolutePath().normalize()
    ARCHIVE_INPUT_MANIFEST(channel.fromPath(manifest_file))
    manifest_rows = channel.fromPath(manifest_file)
        .splitCsv(header: true, sep: '\t')
        .collect()
        .map { rows ->
            if (!rows) error 'manifest has no data rows'
            def root = manifest_file.parent
            def run_ids = [] as Set
            def assembly_files = [:]
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
                def reads = manifest_path(row.fastq, root, 'fastq', row_number)
                def reference = manifest_path(row.reference_fasta, root, 'reference_fasta', row_number)
                def annotation = manifest_path(row.annotation_gtf, root, 'annotation_gtf', row_number)
                def assembly_key = "${species_id}::${assembly_id}"
                def signature = [reference.toString(), annotation.toString()]
                if (assembly_files.containsKey(assembly_key) && assembly_files[assembly_key] != signature)
                    error "manifest rows for ${assembly_key} must use the same reference_fasta and annotation_gtf"
                assembly_files[assembly_key] = signature
                row._reads = reads
                row._reference = reference
                row._annotation = annotation
                row._normalized = [
                    key: assembly_key,
                    species_id: species_id,
                    assembly_id: assembly_id,
                    run_id: run_id,
                    sample_id: safe_id(row.sample_id ?: run_id, 'sample_id', row_number),
                    technology: technology,
                    minimap2_preset: preset
                ]
            }
            rows.collect { row ->
                def meta = row._normalized
                tuple(meta.key, meta, row._reads, row._reference, row._annotation)
            }
        }

    reference_inputs = manifest_rows.flatMap { records ->
        records.groupBy { record -> record[0] }.collect { key, records_for_assembly ->
            def first = records_for_assembly[0]
            tuple(key, first[1], first[3])
        }
    }
    BUILD_REFERENCE_INDEX(reference_inputs)

    alignment_inputs = manifest_rows
        .flatMap { records -> records }
        .combine(BUILD_REFERENCE_INDEX.out.index, by: 0)
        .map { key, meta, reads, _reference, _annotation, _assembly_meta, index, reference_fai ->
            tuple(meta, reads, index, reference_fai)
        }
    ALIGN_SORT_INDEX_QC(alignment_inputs)
    annotation_by_run = manifest_rows.flatMap { records ->
        records.collect { key, meta, _reads, _reference, annotation -> tuple(meta.run_id, annotation) }
    }
    stringtie_inputs = ALIGN_SORT_INDEX_QC.out.aligned
        .map { meta, bam, _bai -> tuple(meta.run_id, meta, bam) }
        .combine(annotation_by_run, by: 0)
        .map { run_id, meta, bam, annotation -> tuple(meta, bam, annotation) }
    STRINGTIE3_ASSEMBLE(stringtie_inputs)

    assembly_groups = STRINGTIE3_ASSEMBLE.out.assembled
        .map { meta, gtf -> tuple(meta.key, meta, gtf) }
        .groupTuple(by: 0)
        .map { key, metas, gtfs -> tuple(key, metas[0], gtfs) }
    assembly_annotations = manifest_rows.flatMap { records ->
        records.groupBy { record -> record[0] }.collect { key, records_for_assembly ->
            tuple(key, records_for_assembly[0][4])
        }
    }
    merge_inputs = assembly_groups
        .combine(assembly_annotations, by: 0)
        .map { key, meta, gtfs, annotation -> tuple(meta, gtfs, annotation) }
    STRINGTIE3_MERGE(merge_inputs)

    quantified_inputs = ALIGN_SORT_INDEX_QC.out.aligned
        .map { meta, bam, _bai -> tuple(meta.key, meta, bam) }
        .combine(STRINGTIE3_MERGE.out.merged.map { meta, merged_gtf -> tuple(meta.key, merged_gtf) }, by: 0)
        .map { key, meta, bam, merged_gtf -> tuple(meta, bam, merged_gtf) }
    STRINGTIE3_QUANTIFY(quantified_inputs)

}
