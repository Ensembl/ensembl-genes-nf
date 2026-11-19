def shell_quote(value: String) {
    "'${value.replace("'", "'\\''")}'"
}

process STAGE_FTP_FILES {
    tag "${meta.assembly_accession}:${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/curl:7.80.0@sha256:a7e3e10f5b162ac7fdf7fd582056037883e8019611d75a138a7808412161e225'
    containerOptions '--bind /etc:/etc'

    input:
    tuple val(meta), path(local_files, stageAs: 'input??')
    val stage_dir

    output:
    tuple val(meta), path('staged'), emit: staged
    path 'input_checksums.tsv', emit: checksums

    script:
    def entries = []
    def local_index = 0
    local_index = local_index + 1
    entries << [kind: 'genome', spec: meta.files.genome, local: true, source: "input${String.format('%02d', local_index)}"]
    meta.files.transcripts.eachWithIndex { spec, i ->
        local_index = local_index + 1
        entries << [kind: "transcript_${String.format('%02d', i + 1)}", spec: spec, local: true, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.stringtie.each { spec ->
        local_index = local_index + 1
        entries << [kind: 'stringtie', spec: spec, local: true, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.scallop.each { spec ->
        local_index = local_index + 1
        entries << [kind: 'scallop', spec: spec, local: true, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.protein_alignments.eachWithIndex { spec, i ->
        local_index = local_index + 1
        entries << [kind: "protein_alignment_${String.format('%02d', i + 1)}", spec: spec, local: true, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.trnascan.each { spec ->
        local_index = local_index + 1
        entries << [kind: 'trnascan', spec: spec, local: true, source: "input${String.format('%02d', local_index)}"]
    }
    meta.files.crams.each { spec -> entries << [kind: spec.name.replace('.cram', ''), spec: spec, local: false] }
    def cram_directory = meta.files.cram_dir ?: ''
    def commands = entries
        .collect { entry ->
            def destination = entry.kind == 'genome' ? "staged/${meta.genome_name}" : entry.kind.startsWith('protein_alignment_') ? "staged/${entry.kind}_${entry.spec.name.replaceFirst(/\.[^.]+$/, '')}.gtf" : "staged/${entry.kind}${entry.kind.startsWith('cram_') ? '.cram' : entry.kind.startsWith('transcript_') ? '.fa' : entry.kind.startsWith('protein_') ? '.fa' : entry.kind == 'stringtie' ? '.gtf' : entry.kind == 'scallop' ? '.gtf' : '.gff'}"
            def checksum = entry.spec.checksum ? "echo ${shell_quote(entry.spec.checksum + '  ' + destination)} | sha256sum -c -" : ':'
            def copy = entry.local && entry.kind == 'genome' && entry.spec.name.toLowerCase().endsWith('.gz') ? "gzip -cd ${shell_quote(entry.source)} > ${shell_quote(destination)}" : entry.local ? "cp -p ${shell_quote(entry.source)} ${shell_quote(destination)}" : "curl -fL --retry 3 --connect-timeout 30 -o ${shell_quote(destination)} -- ${shell_quote(entry.spec.url)}"
            "${copy}\n${checksum}"
        }
        .join('\n')
    def directory_commands = cram_directory
        ? """
cram_names=\$(mktemp)
curl -fL --retry 3 --connect-timeout 30 -- ${shell_quote(cram_directory)} | sed -n 's/.*href=\"\\([^\"]*\\.Aligned\\.SortedByCoordinate\\.out\\.cram\\)\".*/\\1/p' | sort > \"\$cram_names\"
[ -s \"\$cram_names\" ]
cram_index=1
while IFS= read -r cram_name; do
    cram_file=\"staged/cram_\$(printf '%02d' \"\$cram_index\").cram\"
    curl -fL --retry 3 --connect-timeout 30 -o \"\$cram_file\" -- \"${cram_directory}\$cram_name\"
    cram_index=\$((cram_index + 1))
done < \"\$cram_names\"
rm -f \"\$cram_names\"
"""
        : ''
    """
    mkdir -p staged
    ${commands}
    ${directory_commands}
    mkdir -p ${shell_quote("${stage_dir}/${meta.safe_id}")}
    cp -p staged/* ${shell_quote("${stage_dir}/${meta.safe_id}/")}
    sha256sum staged/* > input_checksums.tsv
    """

    stub:
    """
    mkdir -p staged
    touch staged/genome.fa staged/cram_01.cram
    mkdir -p '${stage_dir}/${meta.safe_id}'
    cp -p staged/* '${stage_dir}/${meta.safe_id}/'
    sha256sum staged/* > input_checksums.tsv
    """
}
