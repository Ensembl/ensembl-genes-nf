
process STANDARDISE_BED12 {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pybedtools_pysam_pip_biopython:1dbd8151223e518c"
    
    tag "${meta.id}"

    publishDir "${params.outdir}/standardised_bed12s/${tool}", mode: 'copy'

    input:
    tuple val(meta), val(tool), path(bed_file)
    path genome_fasta
    path genome_fasta_fai

    output:
    tuple val(meta), val(tool), path("*.valid.bed12")

    script:    
    """
    standardise_bed12.py \\
        -i ${bed_file} \\
        -f ${genome_fasta} \\
        --output_prefix ${meta.id} --verbose
    """
    
    stub:
    """
    mkdir -p standardised_bed12s
    touch standardised_bed12s/${meta.id}_0based.bed12
    """
}