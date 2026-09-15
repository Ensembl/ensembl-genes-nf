process VALIDATE_LONG_READ_MANIFEST {
    tag "manifest"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11--he2b4eab_0'

    input:
    path manifest
    path parser

    output:
    path 'normalised_manifest.tsv', emit: manifest
    path 'manifest_validation.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${parser} ${manifest} normalised_manifest.tsv manifest_validation.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    printf 'run_accession\\ttissue\\tfilename\\turl\\tmd5\\tsource\\tplatform\\tdescription\\nSRR000001\\tstub tissue\\tSRR000001.fastq.gz\\thttps://example.org/SRR000001.fastq.gz\\td41d8cd98f00b204e9800998ecf8427e\\tENA\\tPACBIO_SMRT\\tstub\\n' > normalised_manifest.tsv
    printf 'status\\tdetail\\nok\\tstub\\n' > manifest_validation.tsv
    printf '"%s":\\n    python: 3.11.0\\n' '${task.process}' > versions.yml
    """
}
