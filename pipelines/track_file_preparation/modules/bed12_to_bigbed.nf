process BED12_TO_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:482--hdc0a859_0'

    input:
    tuple val(meta), path(bed12), path(chrom_sizes), path(gtf)

    output:
    tuple val(meta), path("${meta.safe_id}.bb"), path(gtf), emit: bigbed

    script:
    """
    bedToBigBed -type=bed12 -tab ${bed12} ${chrom_sizes} ${meta.safe_id}.bb
    test -s ${meta.safe_id}.bb || { echo 'Gene-model BigBed is empty' >&2; exit 1; }
    """

    stub:
    """
    printf 'stub\\n' > ${meta.safe_id}.bb
    """
}
