process BOWTIE_RRNA_FILTER {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::bowtie=1.3.1 bioconda::samtools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f5ca09fd5aab931d9b87c532c69e0122ce5ff8ec88732f906e12108d48425e9/data' :
        'community.wave.seqera.io/library/bowtie_htslib_samtools:e1e242368ffcb5d3' }"

    publishDir "${params.outdir}/rrna_filter", mode: 'copy', pattern: "*_rrna_filter.log"

    input:
    tuple val(meta), path(reads)
    path index  // Bowtie index directory or files

    output:
    tuple val(meta), path("*_no_rrna.fastq.gz"), emit: filtered_fastq
    tuple val(meta), path("*_rrna_filter.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unzip_cmd = reads.name.endsWith('.gz') ? 'zcat' : 'cat'

    """
    # Find the index base name from the staged files
    INDEX=\$(find -L ./ -name "*.1.ebwt" | sed 's/\\.1\\.ebwt\$//')

    # Count input reads
    INPUT_READS=\$(${unzip_cmd} ${reads} | wc -l | awk '{print \$1/4}')

    # Align to rRNA index and extract unmapped reads
    ${unzip_cmd} ${reads} | \\
        bowtie \\
        -p ${task.cpus} \\
        -v 2 \\
        -k 1 \\
        --un ${prefix}_no_rrna.fastq \\
        ${args} \\
        \$INDEX \\
        - \\
        > /dev/null 2> ${prefix}_bowtie.stderr

    # Compress unmapped reads
    gzip ${prefix}_no_rrna.fastq

    # Count output reads (non-rRNA)
    OUTPUT_READS=\$(zcat ${prefix}_no_rrna.fastq.gz | wc -l | awk '{print \$1/4}')

    # Calculate rRNA contamination
    RRNA_READS=\$((\$INPUT_READS - \$OUTPUT_READS))
    PCT_RRNA=\$(awk -v rrna=\$RRNA_READS -v input=\$INPUT_READS 'BEGIN {printf "%.6f", (rrna/input)*100}')

    # Create log file
    cat > ${prefix}_rrna_filter.log <<EOF
Sample: ${prefix}
Input reads: \$INPUT_READS
rRNA reads: \$RRNA_READS
Percentage rRNA: \$PCT_RRNA
Filtered reads (no rRNA): \$OUTPUT_READS
Percentage retained: \$(awk -v out=\$OUTPUT_READS -v input=\$INPUT_READS 'BEGIN {printf "%.6f", (out/input)*100}')
EOF

    # Append bowtie alignment summary
    echo "" >> ${prefix}_rrna_filter.log
    echo "Bowtie alignment summary:" >> ${prefix}_rrna_filter.log
    cat ${prefix}_bowtie.stderr >> ${prefix}_rrna_filter.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(bowtie --version 2>&1 | head -n1 | sed 's/.*version //g')
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_no_rrna.fastq.gz
    cat > ${prefix}_rrna_filter.log <<EOF
Sample: ${prefix}
Input reads: 1000000
rRNA reads: 500000
Percentage rRNA: 50.0
Filtered reads (no rRNA): 500000
Percentage retained: 50.0
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: 1.3.1
        samtools: 1.19
    END_VERSIONS
    """
}
