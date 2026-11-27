process ASSEMBLE_STAGED_FILES {
    tag "${meta.assembly_accession}:${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/curl:7.80.0@sha256:a7e3e10f5b162ac7fdf7fd582056037883e8019611d75a138a7808412161e225'
    containerOptions '--bind /etc:/etc'

    input:
    tuple val(meta), path(local_stage), path(crams)
    val stage_dir

    output:
    tuple val(meta), path('staged'), emit: staged
    path 'input_checksums.tsv', emit: checksums

    script:
    """
    mkdir -p staged
    cp -p '${local_stage}'/* staged/
    cp -p cram_*.cram staged/
    mkdir -p '${stage_dir}/${meta.safe_id}'
    cp -p staged/* '${stage_dir}/${meta.safe_id}/'
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
