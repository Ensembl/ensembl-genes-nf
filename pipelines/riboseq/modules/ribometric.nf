process RIBOMETRIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::biopython bioconda::pysam"
    container "ghcr.io/jackcurragh/ribometric:latest"

    publishDir "${params.outdir}/RiboMetric", mode: 'copy', pattern: "*RiboMetric.{html,json,csv}"
    publishDir "${params.outdir}/RiboMetric/offsets", mode: 'copy', pattern: "*.offsets.tsv"

    errorStrategy 'ignore'

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
    path ribometric_annotation
    path offset_file  // Optional: external offset file (e.g., calculated offsets)

    output:
    tuple val(meta), path("*RiboMetric.html"), emit: html
    tuple val(meta), path("*RiboMetric.json"), emit: json
    tuple val(meta), path("*RiboMetric.csv"),  emit: csv
    tuple val(meta), path("*.offsets.tsv"), optional: true, emit: offsets
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_size = params.ribometric_sample_size ?: 10000000

    // Determine offset configuration
    // Priority: 1) external offset file, 2) global offset, 3) calculation method
    def offset_args = ''
    if (offset_file && offset_file.name != 'NO_OFFSET_FILE') {
        offset_args = "--offset-read-length ${offset_file}"
    } else if (params.ribometric_offset_global) {
        offset_args = "--offset-global ${params.ribometric_offset_global}"
    } else {
        def offset_method = params.ribometric_offset_method ?: 'tripsviz'
        offset_args = "--offset-calculation-method ${offset_method}"
    }
    // Only output offsets when calculating internally (not when using an external file)
    def output_offsets_arg = (offset_file && offset_file.name != 'NO_OFFSET_FILE') ? '' : "--output-offsets ${prefix}.offsets.tsv"
    """
    RiboMetric run \\
        --bam ${transcriptome_bam} \\
        --annotation ${ribometric_annotation} \\
        --threads $task.cpus \\
        --html \\
        --json \\
        --csv \\
        ${offset_args} \\
        ${output_offsets_arg} \\
        -S ${sample_size} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: \$(RiboMetric --version 2>&1 | sed 's/RiboMetric version //g' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_RiboMetric.html
    touch ${prefix}_RiboMetric.json
    touch ${prefix}_RiboMetric.csv
    printf 'read_len\toffset\n28\t12\n29\t12\n30\t12\n' > ${prefix}.offsets.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: 1.0.0
    END_VERSIONS
    """
}
