process DOWNLOAD_FASTQ {
    label 'fastq_download'
    tag "${meta.species_id}/${meta.assembly_id}/${meta.run_id}"
    container 'quay.io/biocontainers/fastq-dl@sha256:9a3fc95f1a287912ac33f678156ae4cb97243cbc8f3f3419521107f0a81fd0a2'

    input:
    tuple val(meta), val(source_uri)

    output:
    tuple val(meta), path('reads.fastq.gz'), emit: reads
    path 'fastq_download.versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    fastq-dl \\
      --accession '${meta.run_id}' \\
      --provider ena \\
      --protocol https \\
      --max-attempts 3 \\
      --cpus ${task.cpus} \\
      --outdir fastq_dl
    mapfile -t fastqs < <(find fastq_dl -maxdepth 1 -type f \\
      \\( -name '*.fastq.gz' -o -name '*.fq.gz' \\) -print | sort)
    if (( \${#fastqs[@]} != 1 )); then
      printf 'expected one FASTQ for %s, found %s\\n' '${meta.run_id}' "\${#fastqs[@]}" >&2
      printf '%s\\n' "\${fastqs[@]}" >&2
      exit 1
    fi
    mv "\${fastqs[0]}" reads.fastq.gz
    test -s reads.fastq.gz
    cat > fastq_download.versions.yml <<VERSIONS
"${task.process}":
    fastq-dl: \$(fastq-dl --version | sed 's/fastq-dl //')
VERSIONS
    """

    stub:
    """
    printf '@stub\\nACGT\\n+\\n!!!!\\n' | gzip -c > reads.fastq.gz
    printf '"DOWNLOAD_FASTQ":\\n    fastq-dl: stub\\n' > fastq_download.versions.yml
    """
}
