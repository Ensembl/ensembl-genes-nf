process VALIDATE_GTF {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(gtf), path(chrom_sizes)

    output:
    tuple val(meta), path(gtf), path(chrom_sizes), emit: validated
    path "${meta.safe_id}.gtf_validation.tsv", emit: validation

    script:
    """
    awk 'NF != 9 { print "Malformed GTF row: " \$0 > "/dev/stderr"; bad=1 } END { if (bad || NR == 0) exit 1 }' ${gtf}
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tgene_model\\tcomplete\\n' > ${meta.safe_id}.gtf_validation.tsv
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tgene_model\\tcomplete\\n' > ${meta.safe_id}.gtf_validation.tsv
    """
}
