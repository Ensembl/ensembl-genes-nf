process STAR_JUNCTIONS_TO_BED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(sj), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.junctions.bed"), path(bam), path(chrom_sizes), emit: bed

    script:
    """
    awk 'BEGIN { OFS="\\t" } NR == 1 && \$1 == "chrom" && \$2 == "start" && \$3 == "end" && \$4 == "strand" && \$5 == "motif" && \$6 == "annotated" && \$7 == "uniqReads" && \$8 == "multiReads" && \$9 == "maxOverhang" { next } { strand=(\$4 == 1 ? "+" : (\$4 == 2 ? "-" : ".")); print \$1, \$2-1, \$3, \$1":"\$2"-"\$3, \$7, strand, \$8 }' ${sj} > ${meta.safe_id}.junctions.bed
    """

    stub:
    """
    printf 'chr1\\t0\\t100\\tchr1:1-100\\t10\\t+\\t0\\n' > ${meta.safe_id}.junctions.bed
    """
}
