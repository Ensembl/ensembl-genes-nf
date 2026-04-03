// DOWNLOAD_ASSEMBLY
// Download genome FASTA and assembly report from NCBI FTP.
// GCA accession (GenBank) path:
//   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/857/705/GCA_001857705.1_Name/

process DOWNLOAD_ASSEMBLY {
    tag "${accession}"
    label 'process_low'

    conda "conda-forge::wget=1.21.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.4--h2b5d971_0' :
        'biocontainers/wget:1.21.4--h2b5d971_0' }"

    input:
    val accession    // GCA_... accession (GenBank)
    val assembly     // assembly name, e.g. GRCh38.p14

    output:
    path "*_genomic.fna.gz",         emit: fasta_gz
    path "*_assembly_report.txt",    emit: report
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = "${accession}_${assembly}"
    def acc_type = accession.substring(0, 3)   // GCA or GCF
    def digits   = accession.replaceAll(/^GC[AF]_/, '').replaceAll(/\.\d+$/, '')
    def p1       = digits[0..2]
    def p2       = digits[3..5]
    def p3       = digits[6..8]
    def base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/${acc_type}/${p1}/${p2}/${p3}/${prefix}"
    """
    wget -q --tries=3 --timeout=120 \\
        "${base_url}/${prefix}_genomic.fna.gz" \\
        -O ${prefix}_genomic.fna.gz

    wget -q --tries=3 --timeout=60 \\
        "${base_url}/${prefix}_assembly_report.txt" \\
        -O ${prefix}_assembly_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //' | cut -d' ' -f1)
    END_VERSIONS
    """

    stub:
    def prefix = "${accession}_${assembly}"
    """
    printf '>chr1\\nATCGATCG\\n' | gzip > ${prefix}_genomic.fna.gz
    printf '# Sequence-Name\\tSequence-Role\\tAssigned-Molecule\\tAssigned-Molecule-Location/Type\\tGenBank-Accn\\tRelationship\\tRefSeq-Accn\\tAssembly-Unit\\tSequence-Length\\tUCSC-style-name\\n' > ${prefix}_assembly_report.txt
    printf 'chr1\\tassembled-molecule\\t1\\tChromosome\\tCM000663.2\\t=\\tNC_000001.11\\tPrimary Assembly\\t248956422\\tchr1\\n' >> ${prefix}_assembly_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.21.4
    END_VERSIONS
    """
}
