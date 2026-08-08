process BAM_TO_BED {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.9 bioconda::pysam=0.23.3 bioconda::samtools=1.20 conda-forge::biopython conda-forge::numpy"
    container "ghcr.io/lapti-ucc/riboseqorg-nf-bam-to-bed:latest"

    input:
    tuple val(meta), path(bam), path(bai), path(offsets)

    output:
    tuple val(meta), path("*.sorted.bedgraph"), emit: bedgraph
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Use BAM filename (without .bam extension) as prefix to preserve filtering type names
    def tier_suffix = meta.track_tier ? ".${meta.track_tier}" : ""
    def prefix = task.ext.prefix ?: "${bam.baseName}${tier_suffix}"
    def is_stranded = args.contains('--stranded')
    """
    bam_to_bed.py \\
        --bam ${bam} \\
        --offsets ${offsets} \\
        --prefix ${prefix} \\
        $args

    if [ "${is_stranded}" == "true" ]; then
        sort -k1,1 -k2,2n ${prefix}.forward.bedgraph > ${prefix}.forward.sorted.bedgraph
        sort -k1,1 -k2,2n ${prefix}.reverse.bedgraph > ${prefix}.reverse.sorted.bedgraph
    else
        sort -k1,1 -k2,2n ${prefix}.bedgraph > ${prefix}.sorted.bedgraph
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        sort: \$(sort --version | sed -n 1p | sed 's/^.*sort (GNU coreutils) //g')
    END_VERSIONS
    """

    stub:
    // Use BAM filename (without .bam extension) as prefix to preserve filtering type names
    def tier_suffix = meta.track_tier ? ".${meta.track_tier}" : ""
    def prefix = task.ext.prefix ?: "${bam.baseName}${tier_suffix}"
    def is_stranded = task.ext.args?.contains('--stranded') ?: false
    """
    if [ "${is_stranded}" == "true" ]; then
        touch ${prefix}.forward.sorted.bedgraph
        touch ${prefix}.reverse.sorted.bedgraph
    else
        touch ${prefix}.sorted.bedgraph
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.5
        pysam: 0.23.3
        sort: 8.32
    END_VERSIONS
    """
}
