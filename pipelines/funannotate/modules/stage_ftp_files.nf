def shell_quote(String value) {
    "'${value.replace("'", "'\\''")}'"
}

process STAGE_FTP_FILES {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    input:
    val meta
    val stage_dir

    output:
    tuple val(meta), path('staged'), emit: staged
    path 'normalized_manifest.tsv', emit: manifest
    path 'input_checksums.tsv', emit: checksums
    path 'acquisition.status.tsv', emit: status

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
        "curl -fL --retry 3 --connect-timeout 30 -- ${shell_quote(entry.spec.url)} -o ${shell_quote(destination)}\n${checksum}\nprintf '%s\\t%s\\t%s\\n' ${shell_quote(entry.kind)} ${shell_quote(destination)} ${shell_quote(entry.spec.url)} >> staged/inputs.tsv"
    }.join('\n')
    def normalized = "cohort_id\\tassembly_accession\\tassembly_release\\tspecies\\tbam_count\\n${meta.id}\\t${meta.assembly_accession}\\t${meta.assembly_release}\\t${meta.species}\\t${meta.bam_count}\\n"
    """
    set -euo pipefail
    mkdir -p staged
    ${commands}
    for input in staged/*; do
        [ "\$input" = staged/inputs.tsv ] || test -s "\$input"
    done
    mkdir -p ${shell_quote("${stage_dir}/${meta.safe_id}")}
    cp -p staged/* ${shell_quote("${stage_dir}/${meta.safe_id}/")}
    printf '%s' '${normalized}' > normalized_manifest.tsv
    sha256sum staged/* > input_checksums.tsv
    printf 'cohort_id\\tstatus\\n%s\\tSTAGED\\n' '${meta.id}' > acquisition.status.tsv
    """

    stub:
    """
    set -euo pipefail
    mkdir -p staged
    touch staged/genome.fa staged/bam_01.bam staged/bam_02.bam
    printf 'genome\\tstaged/genome.fa\\tstub\\n bam_01\\tstaged/bam_01.bam\\tstub\\n bam_02\\tstaged/bam_02.bam\\tstub\\n' > staged/inputs.tsv
    mkdir -p '${stage_dir}/${meta.safe_id}'
    cp -p staged/* '${stage_dir}/${meta.safe_id}/'
    printf 'cohort_id\\tassembly_accession\\tassembly_release\\tspecies\\tbam_count\\n%s\\t%s\\t%s\\t%s\\t%s\\n' '${meta.id}' '${meta.assembly_accession}' '${meta.assembly_release}' '${meta.species}' '${meta.bam_count}' > normalized_manifest.tsv
    printf 'stub\\tSTUB\\n' > input_checksums.tsv
    printf 'cohort_id\\tstatus\\n%s\\tSTAGED\\n' '${meta.id}' > acquisition.status.tsv
    """
}
