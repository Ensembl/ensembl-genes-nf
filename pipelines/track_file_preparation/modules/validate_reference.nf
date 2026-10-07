process VALIDATE_REFERENCE {
    tag { "${assembly_release}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    path chrom_sizes
    val assembly_release

    output:
    path 'reference_check.tsv', emit: checked
    path 'reference.versions.yml', emit: versions

    script:
    """
    echo -e 'assembly_release\tchrom_sizes\tstatus' > reference_check.tsv
    awk 'NF != 2 || \$2 !~ /^[0-9]+\$/ { print "Invalid chrom.sizes row: " \$0 > "/dev/stderr"; bad=1 } END { if (bad || NR == 0) exit 1 }' ${chrom_sizes}
    echo -e '${assembly_release}\t${chrom_sizes}\tcomplete' >> reference_check.tsv
    cat <<-END_VERSIONS > reference.versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    printf 'assembly_release\tchrom_sizes\tstatus\n${assembly_release}\t${chrom_sizes}\tcomplete\n' > reference_check.tsv
    printf '"%s":\n    python: stub\n' '${task.process}' > reference.versions.yml
    """
}
