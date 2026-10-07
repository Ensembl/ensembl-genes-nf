process JUNCTION_BED_TO_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bed), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.bb"), path(bam), emit: bigbed

    script:
    """
    bedToBigBed -type=bed6+1 -tab ${bed} ${chrom_sizes} ${meta.safe_id}.bb
    test -s ${meta.safe_id}.bb || { echo 'Splice-junction BigBed is empty' >&2; exit 1; }
    """

    stub:
    """
    printf 'stub\\n' > ${meta.safe_id}.bb
    """
}
