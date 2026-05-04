/*
 * Merge all study matrices into global Zarr matrix
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
    path "global_matrix.zarr", emit: matrix, type: 'dir'
    path "global_reads.fasta", emit: fasta
    path "global_metadata.parquet", emit: metadata
    path "global_config.json", emit: config
    path "samples.json", emit: samples, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def chunk_size = task.ext.chunk_size ?: 10000
    def partition_flag = task.ext.partition ? "--partition ${task.ext.partition}" : ''
    """
    # Organize flat files into study directories
    # Files are named: {study_id}_matrix.npz, {study_id}_vocab.pkl, {study_id}_sequences.txt.gz
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
        --chunk-size ${chunk_size} \\
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
    mkdir -p global_matrix.zarr
    touch global_reads.fasta
    touch global_metadata.parquet
    echo '{"version": "1.0", "n_reads": 0, "n_samples": 0}' > global_config.json

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
