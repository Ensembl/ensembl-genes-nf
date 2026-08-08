process COLLATE_VERSIONS {
    tag 'software-versions'
    label 'process_light'

    container 'community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2'

    input:
    val version_files

    output:
    path 'software_versions.yml', emit: report

    script:
    def version_args = version_files.collect { value -> "'${value}'" }.join(' ')
    """
    collate_versions.py ${version_args} > software_versions.yml
    """

    stub:
    """
    printf '# No software versions were emitted in stub mode\\n' > software_versions.yml
    """
}
