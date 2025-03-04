process RAPID {
    tag "${meta.id}"
    errorStrategy 'retry'
    maxRetries 0
    publishDir params.fasta_dir, mode: 'copy'

    input:
    tuple val(meta), val(species), val(accession)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: fasta

    script:
    """
    rapid_fetch.py -s '${species}' -a '${accession}' -o '.'
    if [ -f "${meta.id}.fa.gz" ]; then
        gzip -d "${meta.id}.fa.gz"
    fi
    """
}