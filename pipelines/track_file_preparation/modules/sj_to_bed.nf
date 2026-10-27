process STAR_JUNCTIONS_TO_BED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(sj), path(bam), path(chrom_sizes)
    path converter

    output:
    tuple val(meta), path("${meta.safe_id}.junctions.bed"), path(bam), path(chrom_sizes), emit: bed

    script:
    """
    python3 "${converter}" "${sj}" ${meta.safe_id}.junctions.bed
    """

    stub:
    """
    printf 'chr1\\t0\\t100\\tchr1:1-100:+\\t10\\t+\\t10\\t0\\tGT-AG\\n' > ${meta.safe_id}.junctions.bed
    """
}
