def shell_quote(String value) {
    "'${value.replace("'", "'\\''")}'"
}

process STAGE_FTP_FILES {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container params.stage_container ?: 'docker.io/curlimages/curl:8.12.1@sha256:88a9abad9d958340e48564f9bdcdaa29916a2984be59314da709f8bbc0eef6f7'
    input:
    val meta
    val stage_dir

    output:
    tuple val(meta), path('staged'), emit: staged
    path 'input_checksums.tsv', emit: checksums

    script:
    def entries = []
    entries << [kind: 'genome', spec: meta.files.genome]
    meta.files.bams.eachWithIndex { spec, i -> entries << [kind: "bam_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.bais.eachWithIndex { spec, i -> entries << [kind: "bai_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.transcripts.eachWithIndex { spec, i -> entries << [kind: "transcript_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.proteins.eachWithIndex { spec, i -> entries << [kind: "protein_${String.format('%02d', i + 1)}", spec: spec] }
    meta.files.pasa_gff.each { entries << [kind: 'pasa', spec: it] }
    meta.files.parameters_json.each { entries << [kind: 'parameters', spec: it] }
    def commands = entries.collect { entry ->
        def destination = "staged/${entry.kind}${entry.kind == 'genome' ? (entry.spec.name.endsWith('.gz') ? '.fa.gz' : '.fa') : entry.kind.startsWith('bam_') ? '.bam' : entry.kind.startsWith('bai_') ? '.bam.bai' : entry.kind.startsWith('transcript_') ? '.fa' : entry.kind.startsWith('protein_') ? '.fa' : entry.kind == 'pasa' ? '.gff3' : '.json'}"
        def checksum = entry.spec.checksum ? "printf '%s  %s\\n' ${shell_quote(entry.spec.checksum)} ${shell_quote(destination)} | sha256sum -c -" : ':'
        "curl -fL --retry 3 --connect-timeout 30 -- ${shell_quote(entry.spec.url)} -o ${shell_quote(destination)}\n${checksum}"
    }.join('\n')
    """
    set -euo pipefail
    mkdir -p staged
    ${commands}
    mkdir -p ${shell_quote("${stage_dir}/${meta.safe_id}")}
    cp -p staged/* ${shell_quote("${stage_dir}/${meta.safe_id}/")}
    sha256sum staged/* > input_checksums.tsv
    """

    stub:
    """
    set -euo pipefail
    mkdir -p staged
    touch staged/genome.fa staged/bam_01.bam staged/bam_02.bam
    mkdir -p '${stage_dir}/${meta.safe_id}'
    cp -p staged/* '${stage_dir}/${meta.safe_id}/'
    sha256sum staged/* > input_checksums.tsv
    """
}
