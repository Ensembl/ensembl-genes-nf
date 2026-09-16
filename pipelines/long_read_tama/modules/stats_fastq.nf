process STATS_FASTQ {
    tag "${meta.id}"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(reads)
    path inspector

    output:
    tuple val(meta), path('read_validation.tsv'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    test -f "${reads}" || { echo "Expected FASTQ file for ${meta.id}" >&2; exit 1; }
    gzip -t "${reads}" || { echo "Invalid gzip FASTQ for ${meta.id}: ${reads}" >&2; exit 1; }
    python3 "${inspector}" stats "${reads}" read_validation.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    printf 'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\nreads.fastq.gz\tFASTQ\tDNA\t1\t4\t4\t4\t4\n' > read_validation.tsv
    printf '"%s":\n    python: 3.11.0\n' '${task.process}' > versions.yml
    """
}
