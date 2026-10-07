process SORT_JUNCTION_BED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bed), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.sorted.junctions.bed"), path(bam), path(chrom_sizes), emit: sorted

    script:
    """
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n ${bed} > ${meta.safe_id}.sorted.junctions.bed
    """

    stub:
    """
    cp ${bed} ${meta.safe_id}.sorted.junctions.bed
    """
}
