process BUILD_COVERAGE_MATRIX {
    tag { "${meta.gca_accession}:${meta.coverage_strand}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0'
    publishDir "${params.outdir}/coverage", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bigwigs), val(labels)

    output:
    path "${meta.safe_id}.coverage_${meta.coverage_strand}.npz", emit: npz
    path "${meta.safe_id}.coverage_${meta.coverage_strand}.tsv", emit: matrix
    path "${meta.safe_id}.coverage_${meta.coverage_strand}.descriptor.tsv", emit: descriptor
    path "${meta.safe_id}.coverage_${meta.coverage_strand}.versions.yml", emit: versions

    script:
    def bw_args = bigwigs.collect { bigwig -> "'${bigwig}'" }.join(' ')
    def label_args = labels.collect { label -> "'${label}'" }.join(' ')
    """
    multiBigwigSummary bins \
        --bwfiles ${bw_args} \
        --labels ${label_args} \
        --binSize ${meta.summary_bin_size} \
        --numberOfProcessors ${task.cpus} \
        --outFileName ${meta.safe_id}.coverage_${meta.coverage_strand}.npz \
        --outRawCounts ${meta.safe_id}.coverage_${meta.coverage_strand}.tsv
    printf 'output_kind\toutput_file\trun_ids\ncoverage_${meta.coverage_strand}\t${meta.safe_id}.coverage_${meta.coverage_strand}.tsv\t%s\ncoverage_${meta.coverage_strand}\t${meta.safe_id}.coverage_${meta.coverage_strand}.npz\t%s\n' '${labels.join(',')}' '${labels.join(',')}' > ${meta.safe_id}.coverage_${meta.coverage_strand}.descriptor.tsv
    printf '"%s":\n    multiBigwigSummary: \$(multiBigwigSummary --version 2>&1 | sed "s/multiBigwigSummary //")\n' '${task.process}' > ${meta.safe_id}.coverage_${meta.coverage_strand}.versions.yml
    """

    stub:
    def zero_values = labels.collect { _label -> '0' }.join('\\t')
    """
    printf 'chr\tstart\tend\t${labels.join('\t')}\nchr1\t0\t${meta.summary_bin_size}\t${zero_values}\n' > ${meta.safe_id}.coverage_${meta.coverage_strand}.tsv
    touch ${meta.safe_id}.coverage_${meta.coverage_strand}.npz
    printf 'output_kind\toutput_file\trun_ids\ncoverage_${meta.coverage_strand}\t${meta.safe_id}.coverage_${meta.coverage_strand}.tsv\t%s\ncoverage_${meta.coverage_strand}\t${meta.safe_id}.coverage_${meta.coverage_strand}.npz\t%s\n' '${labels.join(',')}' '${labels.join(',')}' > ${meta.safe_id}.coverage_${meta.coverage_strand}.descriptor.tsv
    printf '"%s":\n    multiBigwigSummary: stub\n' '${task.process}' > ${meta.safe_id}.coverage_${meta.coverage_strand}.versions.yml
    """
}
