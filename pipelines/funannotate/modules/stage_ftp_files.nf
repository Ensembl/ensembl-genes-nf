def shell_quote(String value) {
    "'${value.replace("'", "'\\''")}'"
}

process STAGE_FTP_FILES {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container params.stage_container ?: 'docker.io/curlimages/curl:8.12.1@sha256:88a9abad9d958340e48564f9bdcdaa29916a2984be59314da709f8bbc0eef6f7'
    input:
    tuple val(meta), path(genome)
    val stage_dir

    output:
    tuple val(meta), path('staged'), emit: staged
    path 'input_checksums.tsv', emit: checksums

    script:
    def entries = []
    meta.files.crams.eachWithIndex { spec, i -> entries << [kind: "cram_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.crais.eachWithIndex { spec, i -> entries << [kind: "crai_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.transcripts.eachWithIndex { spec, i -> entries << [kind: "transcript_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.proteins.eachWithIndex { spec, i -> entries << [kind: "protein_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.pasa_gff.each { entries << [kind: 'pasa', spec: it] }
    meta.files.parameters_json.each { entries << [kind: 'parameters', spec: it] }
    meta.files.stringtie.each { entries << [kind: 'stringtie', spec: it] }
    meta.files.protein_alignments.each { entries << [kind: 'protein_alignments', spec: it] }
    meta.files.trnascan.each { entries << [kind: 'trnascan', spec: it] }
    def commands = entries.collect { entry ->
        def destination = "staged/${entry.kind}${entry.kind.startsWith('cram_') ? '.cram' : entry.kind.startsWith('crai_') ? '.cram.crai' : entry.kind.startsWith('transcript_') ? '.fa' : entry.kind.startsWith('protein_') ? '.fa' : entry.kind == 'pasa' ? '.gff3' : entry.kind == 'parameters' ? '.json' : entry.kind == 'stringtie' ? '.gtf' : entry.kind == 'protein_alignments' ? '.gff3' : '.gff'}"
        def checksum = entry.spec.checksum ? "echo ${shell_quote(entry.spec.checksum + '  ' + destination)} | sha256sum -c -" : ':'
        "curl -fL --retry 3 --connect-timeout 30 -- ${shell_quote(entry.spec.url)} -o ${shell_quote(destination)}\n${checksum}"
    }.join('\n')
    """
    mkdir -p staged
    cp -p "${genome}" staged/${meta.genome_name}
    ${commands}
    mkdir -p ${shell_quote("${stage_dir}/${meta.safe_id}")}
    cp -p staged/* ${shell_quote("${stage_dir}/${meta.safe_id}/")}
    sha256sum staged/* > input_checksums.tsv
    """

    stub:
    """
    mkdir -p staged
    touch staged/${meta.genome_name} staged/cram_01.cram staged/cram_02.cram
    mkdir -p '${stage_dir}/${meta.safe_id}'
    cp -p staged/* '${stage_dir}/${meta.safe_id}/'
    sha256sum staged/* > input_checksums.tsv
    """
}
