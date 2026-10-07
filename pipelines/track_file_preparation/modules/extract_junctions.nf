process EXTRACT_JUNCTIONS {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1'

    input:
    tuple val(meta), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.derived.sj.out.tab"), path(bam), path(chrom_sizes), emit: junctions

    script:
    """
    regtools junctions extract -o ${meta.safe_id}.derived.sj.out.tab ${bam}
    """

    stub:
    """
    printf 'chr1\\t1\\t100\\t1\\t0\\t2\\t10\\t0\\t0\\n' > ${meta.safe_id}.derived.sj.out.tab
    """
}
