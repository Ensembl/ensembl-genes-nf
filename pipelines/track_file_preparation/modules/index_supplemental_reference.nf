process INDEX_SUPPLEMENTAL_REFERENCE {
    tag { 'supplemental-reference' }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    path supplemental_reference, stageAs: 'supplemental_reference.fasta'

    output:
    path 'supplemental_reference.fasta.fai', emit: fai
    path 'supplemental_reference.samtools.versions.yml', emit: versions

    script:
    """
    samtools faidx supplemental_reference.fasta
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > supplemental_reference.samtools.versions.yml
    """

    stub:
    """
    printf 'NC_000000.1\\t1\\t0\\t1\\t2\\n' > supplemental_reference.fasta.fai
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > supplemental_reference.samtools.versions.yml
    """
}
