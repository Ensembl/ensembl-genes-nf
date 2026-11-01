process JUNCTION_BED_TO_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}:splice_junction" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:482--hdc0a859_0'

    input:
    tuple val(meta), path(bed), path(chrom_sizes)
    path autosql

    output:
    tuple val(meta), path("${meta.safe_id}.bb"), emit: bigbed
    path "${meta.safe_id}.junction_bed_to_bigbed.versions.yml", emit: versions

    script:
    """
    bedToBigBed -type=bed6+3 -as=${autosql} -tab ${bed} ${chrom_sizes} ${meta.safe_id}.bb
    test -s ${meta.safe_id}.bb || { echo 'Splice-junction BigBed is empty' >&2; exit 1; }
    printf '"%s":\\n    bedToBigBed: 482\\n' '${task.process}' > ${meta.safe_id}.junction_bed_to_bigbed.versions.yml
    """

    stub:
    """
    printf 'stub\\n' > ${meta.safe_id}.bb
    printf '"%s":\\n    bedToBigBed: stub\\n' '${task.process}' > ${meta.safe_id}.junction_bed_to_bigbed.versions.yml
    """
}
