process BUILD_CROSS_SPECIES {
    tag 'cross-species'
    container 'quay.io/biocontainers/python:3.12.12'
    label 'cage_python'
    publishDir "${params.outdir}/cross_species", mode: 'copy'

    input:
    path atlas
    path orthology_manifest
    path helper

    output:
    path 'cross_species_comparison.tsv', emit: comparison

    script:
    """
    python ${helper} --atlas ${atlas.join(' ')} --orthology ${orthology_manifest} --output cross_species_comparison.tsv
    """

    stub:
    """
    printf 'species\\tassembly\\tortholog_group\\tstatus\\n' > cross_species_comparison.tsv
    """
}
