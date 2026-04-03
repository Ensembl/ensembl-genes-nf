// DOWNLOAD_REFSEQ
// Download the NCBI RefSeq GFF3 for a given assembly accession + name.
// FTP URL pattern:
//   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/857/705/
//     GCF_001857705.1_AssemblyName/GCF_001857705.1_AssemblyName_genomic.gff.gz

process DOWNLOAD_REFSEQ {
    tag "${accession}"
    label 'process_low'

    conda "conda-forge::wget=1.21.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.4--h2b5d971_0' :
        'biocontainers/wget:1.21.4--h2b5d971_0' }"

    input:
    val accession   // e.g. GCF_001857705.1
    val assembly    // e.g. MusMus_1.0

    output:
    path "*_genomic.gff.gz", emit: gff_gz
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Build FTP path from accession: GCF_001857705.1 → GCF/001/857/705
    def acc_digits = accession.replaceAll(/^GCF_/, '').replaceAll(/\.\d+$/, '')
    def p1 = acc_digits[0..2]
    def p2 = acc_digits[3..5]
    def p3 = acc_digits[6..8]
    def p4 = acc_digits[9..-1]
    def base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${p1}/${p2}/${p3}/${accession}_${assembly}"
    def filename = "${accession}_${assembly}_genomic.gff.gz"
    """
    wget -q --tries=3 --timeout=60 \\
        "${base_url}/${filename}" \\
        -O ${filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //' | cut -d' ' -f1)
    END_VERSIONS
    """

    stub:
    def filename = "${accession}_${assembly}_genomic.gff.gz"
    """
    printf '##gff-version 3\\n' | gzip > ${filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.21.4
    END_VERSIONS
    """
}
