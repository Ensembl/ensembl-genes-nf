process ACQUIRE_LONG_READ_BAM {
    tag "${meta.id}:${meta.classification}"
    label 'process_high_memory'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(expected_md5), val(expected_filename), val(source_uri), val(sidecar_uris), val(sidecar_md5s)
    val cache_dir

    output:
    tuple val(meta), path('downloaded/*'), emit: artifacts
    path 'versions.yml', emit: versions

    script:
    """
    acquire_bam.py ${expected_md5} ${expected_filename} '${source_uri}' '${sidecar_uris}' '${sidecar_md5s}' '${cache_dir}' '${meta.classification}' '${meta.id}' downloaded downloaded/checksum.tsv
    printf '"%s":\\n    acquisition: curl-and-md5sum\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p downloaded
    touch downloaded/${expected_filename}
    touch downloaded/placeholder.pbi
    printf 'run_accession\\texpected_md5\\tactual_md5\\tfilename\\n${meta.id}\\tstub\\tstub\\t${expected_filename}\\n' > downloaded/checksum.tsv
    printf '"%s":\\n    acquisition: stub\\n' '${task.process}' > versions.yml
    """
}
