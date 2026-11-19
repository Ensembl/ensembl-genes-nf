process CONVERT_PROTEIN_ALIGNMENTS {
    tag "${meta.assembly_accession}:${meta.id}"
    label 'python'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: converted
    path 'versions.yml', emit: versions

    script:
    def input_gtfs = meta.files.protein_alignments ? 'staged/protein_alignment_*.gtf' : ''
    """
    if [[ -n '${input_gtfs}' ]]; then
        convert_protein_alignments.py ${input_gtfs} --output staged/genblast_protein_alignments.gff3
        [[ -s staged/genblast_protein_alignments.gff3 ]]
    fi
    cat > versions.yml <<VERSIONS
"${task.process}":
    python: \$(python --version 2>&1)
VERSIONS
    """

    stub:
    """
    ${meta.files.protein_alignments ? 'touch staged/genblast_protein_alignments.gff3' : ':'}
    touch versions.yml
    """
}
