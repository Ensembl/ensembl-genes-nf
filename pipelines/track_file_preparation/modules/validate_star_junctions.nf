process VALIDATE_STAR_JUNCTIONS {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(sj), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path(sj), path(bam), path(chrom_sizes), emit: validated
    path "${meta.safe_id}.sj_validation.tsv", emit: validation

    script:
    """
    awk 'NF != 9 || \$2 < 1 || \$3 < \$2 { print "Malformed STAR junction row: " \$0 > "/dev/stderr"; bad=1 } END { if (bad || NR == 0) exit 1 }' ${sj}
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tsplice_junction\\tcomplete\\n' > ${meta.safe_id}.sj_validation.tsv
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tsplice_junction\\tcomplete\\n' > ${meta.safe_id}.sj_validation.tsv
    """
}
