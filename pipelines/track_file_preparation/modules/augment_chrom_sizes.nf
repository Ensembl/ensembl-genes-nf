process AUGMENT_CHROM_SIZES {
    tag { 'supplemental-reference' }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

    input:
    path chrom_sizes
    path supplemental_fai
    path builder

    output:
    path 'augmented.chrom.sizes', emit: chrom_sizes
    path 'supplemental_sequences.tsv', emit: report
    path 'augment_chrom_sizes.versions.yml', emit: versions

    script:
    """
    python3 ${builder} \\
        --chrom-sizes ${chrom_sizes} \\
        --supplemental-fai ${supplemental_fai} \\
        --output augmented.chrom.sizes \\
        --report supplemental_sequences.tsv \\
        --versions augment_chrom_sizes.versions.yml
    """

    stub:
    """
    cp ${chrom_sizes} augmented.chrom.sizes
    printf 'sequence_name\\tlength\\tsource\\tstatus\\n' > supplemental_sequences.tsv
    printf '"%s":\\n    python: stub\\n' '${task.process}' > augment_chrom_sizes.versions.yml
    """
}
