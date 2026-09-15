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
    set -euo pipefail
    mkdir -p downloaded "${cache_dir}/${meta.classification}/${expected_md5}"
    cache_file="${cache_dir}/${meta.classification}/${expected_md5}/${expected_filename}"
    cache_lock="\${cache_file}.lock"
    while ! mkdir "\${cache_lock}" 2>/dev/null; do sleep 10; done
    trap 'rmdir "\${cache_lock}" 2>/dev/null || true' EXIT

    if [ -f "\${cache_file}" ]; then
        gzip -t "\${cache_file}"
        actual=\$(md5sum "\${cache_file}" | awk '{print \$1}')
        test "\${actual}" = "${expected_md5}"
    else
        mkdir -p acquisition
        if [ -n "${source_uri}" ]; then
            curl --fail --location --retry 3 --output "acquisition/${expected_filename}" "${source_uri}"
        else
            (
                cd acquisition
                fastq-dl -a ${run} --cpus ${task.cpus} ${task.ext.args ?: ''}
            )
        fi
        mapfile -t files < <(find acquisition -type f -name '*.fastq.gz' -print)
        if [ "\${#files[@]}" -ne 1 ]; then
            echo "Expected one long-read FASTQ for ${run}, found \${#files[@]}" >&2
            exit 1
        fi
        gzip -t "\${files[0]}"
        test "\$(basename "\${files[0]}")" = "${expected_filename}"
        actual=\$(md5sum "\${files[0]}" | awk '{print \$1}')
        test "\${actual}" = "${expected_md5}"
        mv "\${files[0]}" "\${cache_file}"
    fi
    # Keep a task-local output for Nextflow while the cache remains persistent.
    cp -p "\${cache_file}" "downloaded/${expected_filename}"
    printf 'run_accession\\texpected_md5\\tactual_md5\\tfilename\\n${run}\\t${expected_md5}\\t\${actual}\\t${expected_filename}\\n' > downloaded/checksum.tsv
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
