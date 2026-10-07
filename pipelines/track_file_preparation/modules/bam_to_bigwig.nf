process BAM_TO_BIGWIG {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.bw"), path(bam), emit: bigwig
    path "${meta.safe_id}.coverage.versions.yml", emit: versions

    script:
    def bin_size = task.ext.bin_size ?: params.coverage_bin_size
    def normalization = task.ext.normalization ?: params.coverage_normalization
    """
    bamCoverage --bam ${bam} --outFileName ${meta.safe_id}.bw --binSize ${bin_size} --normalizeUsing ${normalization} --numberOfProcessors ${task.cpus} --effectiveGenomeSize 0
    test -s ${meta.safe_id}.bw || { echo 'Coverage BigWig is empty' >&2; exit 1; }
    printf '"%s":\\n    deepTools: %s\\n' '${task.process}' "\$(bamCoverage --version 2>&1 | head -n1)" > ${meta.safe_id}.coverage.versions.yml
    """

    stub:
    """
    printf 'stub\\n' > ${meta.safe_id}.bw
    printf '"%s":\\n    deepTools: stub\\n' '${task.process}' > ${meta.safe_id}.coverage.versions.yml
    """
}
