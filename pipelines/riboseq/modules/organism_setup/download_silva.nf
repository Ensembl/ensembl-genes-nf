/*
 * DOWNLOAD_SILVA
 * Download rRNA sequences from SILVA database
 */

process DOWNLOAD_SILVA {
    tag "SILVA"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' }"

    publishDir "${params.outdir}/organism_setup/${organism}/${version}", mode: 'copy'

    input:
    val(silva_url)
    val(organism)
    val(version)

    output:
    path "rRNA.fa",        emit: rrna_fasta
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Download SILVA rRNA database
    # Common URLs:
    # - SSU (16S/18S): https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
    # - LSU (23S/28S): https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz

    curl --silent --fail -L -o silva_rrna.fasta.gz "${silva_url}"

    # Decompress
    gzip -d silva_rrna.fasta.gz

    # Rename to standard output name
    mv silva_rrna.fasta rRNA.fa

    echo "Downloaded SILVA rRNA sequences: \$(grep -c '^>' rRNA.fa) sequences"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n1 | sed 's/curl //' | cut -d' ' -f1)
    END_VERSIONS
    """

    stub:
    """
    touch rRNA.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: 8.4.0
    END_VERSIONS
    """
}
