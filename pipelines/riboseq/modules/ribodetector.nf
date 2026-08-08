process RIBODETECTOR {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::ribodetector=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribodetector:0.3.1--pyhdfd78af_0' :
        'biocontainers/ribodetector:0.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_no_rrna.fastq.gz"), emit: filtered_fastq
    tuple val(meta), path("*_rrna.fastq.gz"), emit: rrna_fastq
    tuple val(meta), path("*_rrna_filter.log"), emit: log
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unzip_cmd = reads.name.endsWith('.gz') ? 'zcat' : 'cat'
    def chunk_size = params.ribodetector_chunk_size ?: 256
    def read_len = params.ribodetector_len ?: 100

    """
    # Count input reads
    INPUT_READS=\$(${unzip_cmd} ${reads} | wc -l | awk '{print \$1/4}')

    # Run RiboDetector
    # -i: input file
    # -o: output file (non-rRNA reads)
    # -r: output file for rRNA reads
    # -t: number of threads
    # -l: read length (for short reads, use actual length)
    # -c: chunk size for processing
    ribodetector_cpu \\
        -i ${reads} \\
        -o ${prefix}_no_rrna.fastq.gz \\
        -r ${prefix}_rrna.fastq.gz \\
        -t ${task.cpus} \\
        -l ${read_len} \\
        -c ${chunk_size} \\
        ${args}

    # Count output reads (non-rRNA)
    OUTPUT_READS=\$(zcat ${prefix}_no_rrna.fastq.gz | wc -l | awk '{print \$1/4}')

    # Count rRNA reads
    RRNA_READS=\$(zcat ${prefix}_rrna.fastq.gz | wc -l | awk '{print \$1/4}')

    # Calculate percentages
    PCT_RRNA=\$(awk -v rrna=\$RRNA_READS -v input=\$INPUT_READS 'BEGIN {printf "%.6f", (rrna/input)*100}')

    # Create log file
    cat > ${prefix}_rrna_filter.log <<EOF
Sample: ${prefix}
Method: RiboDetector (ML-based)
Input reads: \$INPUT_READS
rRNA reads: \$RRNA_READS
Percentage rRNA: \$PCT_RRNA
Filtered reads (no rRNA): \$OUTPUT_READS
Percentage retained: \$(awk -v out=\$OUTPUT_READS -v input=\$INPUT_READS 'BEGIN {printf "%.6f", (out/input)*100}')
Parameters:
  Read length: ${read_len}
  Chunk size: ${chunk_size}
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: \$(ribodetector_cpu --version 2>&1 | sed 's/ribodetector_cpu //g' || echo "0.3.1")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_no_rrna.fastq.gz
    touch ${prefix}_rrna.fastq.gz
    cat > ${prefix}_rrna_filter.log <<EOF
Sample: ${prefix}
Method: RiboDetector (ML-based)
Input reads: 1000000
rRNA reads: 450000
Percentage rRNA: 45.0
Filtered reads (no rRNA): 550000
Percentage retained: 55.0
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: 0.3.1
    END_VERSIONS
    """
}
