process MERGE_UNIQUE_READS {
    tag "${mode}"
    label 'process_high'

    conda "conda-forge::python=3.10 conda-forge::polars=0.20.0"
    container 'community.wave.seqera.io/library/python_polars:356ea2c3bf4ab1ed' 

    publishDir "${params.outdir}/unique_reads", mode: 'copy', pattern: "unique_reads*"
    publishDir "${params.outdir}/unique_reads", mode: 'copy', pattern: "count_matrix*"
    publishDir "${params.outdir}/unique_reads", mode: 'copy', pattern: "read_mapping*"
    publishDir "${params.outdir}/unique_reads", mode: 'copy', pattern: "processing_summary*"

    input:
    path collapsed_fastas       // List of all collapsed FASTA files from current run
    path previous_index         // Optional: previous unique_reads.fa, count_matrix.parquet, read_mapping.json
    val mode                    // 'per_run' or 'progressive'

    output:
    path "unique_reads*.fa", emit: unique_reads_fasta
    path "count_matrix*.parquet", emit: count_matrix
    path "read_mapping*.json", emit: read_mapping
    path "processing_summary*.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def batch_size = task.ext.batch_size ?: 20
    def n_workers = task.ext.n_workers ?: task.cpus
    def progressive_args = mode == 'progressive' && previous_index ? "--previous-index ${previous_index}" : ""
    def suffix = mode == 'per_run' ? "_run_${workflow.runName}" : ""
    """
    # Create a file list for the Python script
    find . -name "*.fa" -o -name "*.fa.gz" > file_list.txt

    # Run the unique reads merger
    sorted_read_index_merger.py \\
        --file-list file_list.txt \\
        --output-dir . \\
        --batch-size ${batch_size} \\
        --n-workers ${n_workers} \\
        ${progressive_args} \\
        --output-suffix "${suffix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch unique_reads.fa
    touch count_matrix.parquet
    touch read_mapping.json
    touch processing_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: unknown
        polars: unknown
    END_VERSIONS
    """
}
