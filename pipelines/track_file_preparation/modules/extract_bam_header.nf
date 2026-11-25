process EXTRACT_BAM_HEADER {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.gca_accession}.${meta.safe_id}.header.sam"), path(chrom_sizes), emit: header
    path "${meta.gca_accession}.${meta.safe_id}.samtools.versions.yml", emit: versions

    script:
    """
    samtools view -H ${bam} > ${meta.gca_accession}.${meta.safe_id}.header.sam
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > ${meta.gca_accession}.${meta.safe_id}.samtools.versions.yml
    """

    stub:
    """
    printf '@HD\\tVN:1.6\\n@SQ\\tSN:stub\\tLN:1\\n' > ${meta.gca_accession}.${meta.safe_id}.header.sam
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > ${meta.gca_accession}.${meta.safe_id}.samtools.versions.yml
    """
}
