def shell_quote(value: String) {
    "'${value.replace("'", "'\\''")}'"
}

process STAGE_LOCAL_FILES {
    tag "${meta.assembly_accession}:${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/curl:7.80.0@sha256:a7e3e10f5b162ac7fdf7fd582056037883e8019611d75a138a7808412161e225'
    containerOptions '--bind /etc:/etc'

    input:
    tuple val(meta), path(local_files, stageAs: 'input??')

    output:
    tuple val(meta), path('staged_local'), emit: staged

    script:
    def entries = []
    def local_index = 0
    local_index = local_index + 1
    entries << [kind: 'genome', spec: meta.files.genome, source: "input${String.format('%02d', local_index)}"]
    meta.files.transcripts.eachWithIndex { spec, i ->
        local_index = local_index + 1
        entries << [kind: "transcript_${String.format('%02d', i + 1)}", spec: spec, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.stringtie.each { spec ->
        local_index = local_index + 1
        entries << [kind: 'stringtie', spec: spec, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.scallop.each { spec ->
        local_index = local_index + 1
        entries << [kind: 'scallop', spec: spec, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.protein_alignments.eachWithIndex { spec, i ->
        local_index = local_index + 1
        entries << [kind: "protein_alignment_${String.format('%02d', i + 1)}", spec: spec, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.trnascan.each { spec ->
        local_index = local_index + 1
        entries << [kind: 'trnascan', spec: spec, source: "input${String.format('%02d', local_index)}"]
    }
    def commands = entries
        .collect { entry ->
            def destination = entry.kind == 'genome' ? 'staged_local/genome.fa' : entry.kind.startsWith('protein_alignment_') ? "staged_local/${entry.kind}_${entry.spec.name.replaceFirst(/\.[^.]+$/, '')}.gtf" : "staged_local/${entry.kind}${entry.kind.startsWith('transcript_') ? '.fa' : entry.kind == 'stringtie' || entry.kind == 'scallop' ? '.gtf' : '.gff'}"
            def copy = entry.kind == 'genome' && entry.spec.name.toLowerCase().endsWith('.gz') ? "gzip -cd ${shell_quote(entry.source)} > ${shell_quote(destination)}" : "cp -p ${shell_quote(entry.source)} ${shell_quote(destination)}"
            "${copy}"
        }
        .join('\n')
    """
    mkdir -p staged_local
    ${commands}
    """

    stub:
    """
    mkdir -p staged_local
    touch staged_local/genome.fa
    """
}
