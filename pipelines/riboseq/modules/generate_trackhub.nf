/*
 * GENERATE_TRACKHUB
 * Generate UCSC-compatible track hub from BigWig files
 */

process GENERATE_TRACKHUB {
    tag "${hub_name}"
    label 'process_low'

    container "community.wave.seqera.io/library/pip_trackhub:b1b9686e5cada428"

    publishDir "${params.outdir}/trackhubs", mode: 'copy'

    input:
    path(bigwig_files)
    val(hub_name)
    val(genome)
    val(email)
    val(sample_regex)
    val(annotation_regex)

    output:
    path "trackhub_output/**",   emit: trackhub
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bigwig_paths = bigwig_files ? "--bigwig " + bigwig_files.collect { "'$it'" }.join(' ') : ''
    def sample_regex_param = sample_regex ? "--sample-regex '${sample_regex}'" : ''
    def annotation_regex_param = annotation_regex ? "--annotation-regex '${annotation_regex}'" : ''
    """
    TrackHubGenerator.py create \\
        ${bigwig_paths} \\
        --hub-name "${hub_name}" \\
        --genome "${genome}" \\
        --output-dir trackhub_output \\
        ${sample_regex_param} \\
        ${annotation_regex_param} \\
        --email "${email}" \\
        ${args}

    # Version reporting
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TrackHubGenerator: "1.0.0"
        python: \$(python --version | sed 's/Python //')
        trackhub: \$(python -c "import trackhub; print(trackhub.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p trackhub_output/${hub_name}
    touch trackhub_output/${hub_name}/hub.txt
    touch trackhub_output/${hub_name}/genomes.txt
    mkdir -p trackhub_output/${hub_name}/${genome}
    touch trackhub_output/${hub_name}/${genome}/trackDb.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TrackHubGenerator: "1.0.0"
        python: 3.10
        trackhub: 0.2.0
    END_VERSIONS
    """
}
