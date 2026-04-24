// FETCH_READS_FROM_ENA
// Download short-read RNA-seq data from ENA by BioProject or run accessions.
// Outputs a sample sheet CSV (id,fastq_1[,fastq_2],strandedness) and the
// downloaded FASTQ files, ready for STAR alignment.
//
// Strategy:
//   1. Query the ENA Portal API for run accessions belonging to the BioProject/study
//   2. For each run, query sample metadata (tissue/cell_type for id field)
//   3. Download paired FASTQ from ENA FTP (or fall back to SRA if FTP fails)
//   4. Write sample_sheet.csv
//
// Accession formats accepted:
//   PRJNA123456   (NCBI BioProject)
//   PRJEB123456   (ENA BioProject)
//   SRR123456,SRR123457  (individual run accessions, comma-separated)
//   ERR123456,ERR123457

process FETCH_READS_FROM_ENA {
    tag "${accession}"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::requests=2.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/raw_reads", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.outdir}",           mode: 'copy', pattern: "sample_sheet.csv"

    input:
    val accession   // BioProject or comma-separated run accessions

    output:
    path "sample_sheet.csv",  emit: sample_sheet
    path "*.fastq.gz",        emit: fastq,      optional: true
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_runs = params.max_rnaseq_runs ?: 50   // cap downloads to avoid huge datasets
    def strandedness = params.rnaseq_strandedness ?: 'auto'
    """
    fetch_reads_from_ena.py \\
        --accession   "${accession}" \\
        --outdir      . \\
        --max-runs    ${max_runs} \\
        --strandedness ${strandedness}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """

    stub:
    """
    # Simulate two paired-end samples
    printf 'id,fastq_1,fastq_2,strandedness\\n' > sample_sheet.csv
    printf 'liver,liver_R1.fastq.gz,liver_R2.fastq.gz,unstranded\\n' >> sample_sheet.csv
    printf 'kidney,kidney_R1.fastq.gz,kidney_R2.fastq.gz,unstranded\\n' >> sample_sheet.csv

    touch liver_R1.fastq.gz liver_R2.fastq.gz kidney_R1.fastq.gz kidney_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        requests: 2.31.0
    END_VERSIONS
    """
}
