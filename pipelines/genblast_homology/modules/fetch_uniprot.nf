// FETCH_UNIPROT
// Download reviewed UniProt proteins for a given taxon (clade or species).
// Uses the UniProt REST API to stream sequences in FASTA format.
//
// The taxon_id should be the NCBI taxonomy ID of the clade you want
// (e.g. 9989 for Rodentia, 40674 for Mammalia, 7742 for Vertebrata).
// Reviewed (Swiss-Prot) proteins only are downloaded to keep the set
// high-quality and manageable.
//
// Output: <taxon_id>_reviewed.fa

process FETCH_UNIPROT {
    tag "${taxon_id}"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::requests=2.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/uniprot", mode: 'copy', pattern: "*.fa"

    input:
    val taxon_id   // NCBI taxonomy ID (integer or string)

    output:
    path "*.fa",       emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fetch_uniprot.py \\
        --taxon-id  ${taxon_id} \\
        --out       ${taxon_id}_reviewed.fa \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """

    stub:
    """
    printf '>sp|Q99999|STUB_MOUSE Stub protein OS=Mus musculus OX=10090\\nMPQSTUBSEQUENCE\\n' \\
        > ${taxon_id}_reviewed.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        requests: 2.31.0
    END_VERSIONS
    """
}
