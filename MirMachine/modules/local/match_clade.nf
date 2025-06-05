

process MATCH_CLADE {
    tag "${meta.id}"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fasta)
    path clade_file
    path taxa_db

    output:
    tuple val(meta), path("${meta.id}_clade_match.txt"), emit: matched_clade

    script:
    """
    match_clade.py -s "${meta.species}" -c $clade_file -o ${meta.id}_clade_match.txt
    """
}
