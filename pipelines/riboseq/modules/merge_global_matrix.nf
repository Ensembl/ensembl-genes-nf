/*
 * Merge all study matrices into a global count store
 * Outputs: global matrix, unique reads FASTA for alignment, metadata
 */

process MERGE_GLOBAL_MATRIX {
    tag "global"
    label 'process_ultra_high'

    // Dependencies: python, numpy, scipy, zarr, numcodecs, polars, xxhash
    // Use conda profile or enable Wave for automatic container generation
    conda "conda-forge::python=3.11 conda-forge::numpy=1.26 conda-forge::scipy=1.12 conda-forge::zarr=2.18 conda-forge::numcodecs=0.12 conda-forge::polars=0.20 conda-forge::xxhash-python=3.4"
    container "oras://community.wave.seqera.io/library/pip_numpy_polars_scipy_pruned:ca114eb799eb08b3"
    publishDir "${params.outdir}/global", mode: 'copy'

    input:
    path study_files  // Collected study files (flat: *_matrix.npz, *_vocab.pkl, *_sequences.txt.gz)

    output:
    path "global_matrix.*", emit: matrix
    path "global_counts", emit: counts, optional: true
    path "global_reads.parquet", emit: reads, optional: true
    path "global_samples.parquet", emit: samples_parquet, optional: true
    path "global_studies.parquet", emit: studies_parquet, optional: true
    path "global_tombstones.parquet", emit: tombstones, optional: true
    path "global_retained_counts.parquet", emit: retained_counts, optional: true
    path "global_reads.fasta", emit: fasta
    path "global_metadata.parquet", emit: metadata
    path "global_config.json", emit: config
    path "global_matrix_manifest.json", emit: manifest, optional: true
    path "samples.json", emit: samples, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def chunk_size = task.ext.chunk_size ?: 10000
    def metadata_shard_rows = task.ext.metadata_shard_rows ?: 100000000
    def matrix_format = task.ext.matrix_format ?: 'sparse-parquet'
    def sparse_shard_rows = task.ext.sparse_shard_rows ?: 5000000
    def sparse_read_bucket_size = task.ext.sparse_read_bucket_size ?: 100000
    def append_flag = task.ext.append_to ? "--append-to ${task.ext.append_to}" : ''
    def fasta_flag = task.ext.write_fasta == false ? '--no-write-fasta' : ''
    def partition_flag = task.ext.partition ? "--partition ${task.ext.partition}" : ''
    """
    # Organize flat files into study directories
    # Files are named: {study_id}_matrix.npz, {study_id}_sequences.txt.gz, {study_id}_metadata.json
    for matrix in *_matrix.npz; do
        study_id=\$(basename \$matrix _matrix.npz)
        mkdir -p studies/\$study_id
        mv \${study_id}_matrix.npz studies/\$study_id/ 2>/dev/null || true
        mv \${study_id}_vocab.pkl studies/\$study_id/ 2>/dev/null || true
        mv \${study_id}_sequences.txt.gz studies/\$study_id/ 2>/dev/null || true
        mv \${study_id}_metadata.json studies/\$study_id/ 2>/dev/null || true
    done

    merge_global_matrix.py \\
        --study-dirs studies/* \\
        --output-dir . \\
        --matrix-format ${matrix_format} \\
        --chunk-size ${chunk_size} \\
        --metadata-shard-rows ${metadata_shard_rows} \\
        --sparse-shard-rows ${sparse_shard_rows} \\
        --sparse-read-bucket-size ${sparse_read_bucket_size} \\
        ${append_flag} \\
        ${fasta_flag} \\
        ${partition_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        zarr: \$(python -c "import zarr; print(zarr.__version__)" 2>/dev/null || echo "not installed")
        polars: \$(python -c "import polars; print(polars.__version__)" 2>/dev/null || echo "not installed")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p global_matrix.parquet
    mkdir -p global_counts
    touch global_reads.parquet
    touch global_samples.parquet
    touch global_studies.parquet
    touch global_tombstones.parquet
    touch global_retained_counts.parquet
    touch global_reads.fasta
    touch global_metadata.parquet
    echo '{"version": "1.0", "matrix_format": "sparse-parquet", "n_reads": 0, "n_samples": 0}' > global_config.json
    echo '{"version": "1.0", "matrix_format": "sparse-parquet", "status": "stub"}' > global_matrix_manifest.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        numpy: stub
        scipy: stub
        zarr: stub
        polars: stub
    END_VERSIONS
    """
}
