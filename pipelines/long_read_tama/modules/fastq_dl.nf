process FASTQ_DL {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::fastq-dl=2.0.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-dl:2.0.1--pyhdfd78af_0' :
        'biocontainers/fastq-dl:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(run), val(expected_md5), val(expected_filename), val(source_uri)
    val cache_dir

    output:
    tuple val(meta), path('downloaded/*.fastq.gz'), emit: fastq
    tuple val(meta), path('downloaded/checksum.tsv'), emit: checksum
    path 'versions.yml', emit: versions

    script:
    """
    download_fastq.py ${run} ${expected_md5} ${expected_filename} '${source_uri}' '${cache_dir}' '${meta.classification}' downloaded ${task.cpus} downloaded/checksum.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: \$(fastq-dl --version | sed 's/fastq-dl //')
        md5sum: system
    END_VERSIONS
    """

    stub:
    """
    mkdir -p downloaded
    printf '@stub\\nACGT\\n+\\n!!!!\\n' | gzip -c > downloaded/${expected_filename}
    printf 'run_accession\\texpected_md5\\tactual_md5\\tfilename\\n${run}\\tstub\\tstub\\t${expected_filename}\\n' > downloaded/checksum.tsv
    printf '"%s":\\n    fastq-dl: 2.0.1\\n' '${task.process}' > versions.yml
    """
}
