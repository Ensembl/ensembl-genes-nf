/*
 * URL_DOWNLOAD
 * Download genome and annotation files from custom URLs
 * Supports gzipped files (auto-decompressed)
 */

process URL_DOWNLOAD {
    tag "${organism}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' }"

    publishDir "${params.outdir}/organism_setup/${organism}/${version}", mode: 'copy'

    input:
    val(fasta_url)
    val(gtf_url)
    val(organism)
    val(version)

    output:
    path "*.fa",           emit: genome_fasta
    path "*.gtf",          emit: genome_gtf
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fasta_filename = fasta_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    def gtf_filename = gtf_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    """
    # Download genome FASTA
    curl --silent --fail -L -o genome_download.fa.gz "${fasta_url}" || \\
    curl --silent --fail -L -o genome_download.fa "${fasta_url}"

    # Download GTF annotation
    curl --silent --fail -L -o annotation_download.gtf.gz "${gtf_url}" || \\
    curl --silent --fail -L -o annotation_download.gtf "${gtf_url}"

    # Decompress if gzipped
    if [ -f genome_download.fa.gz ]; then
        gzip -d genome_download.fa.gz
    fi
    if [ -f annotation_download.gtf.gz ]; then
        gzip -d annotation_download.gtf.gz
    fi

    # Rename to clean names
    mv genome_download.fa ${fasta_filename}
    mv annotation_download.gtf ${gtf_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n1 | sed 's/curl //' | cut -d' ' -f1)
    END_VERSIONS
    """

    stub:
    def fasta_filename = fasta_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    def gtf_filename = gtf_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    """
    touch ${fasta_filename}
    touch ${gtf_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: 8.4.0
    END_VERSIONS
    """
}
