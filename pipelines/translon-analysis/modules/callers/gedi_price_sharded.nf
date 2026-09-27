process GEDI_PRICE_SHARDED {
    tag "${meta.id}:${contig}"
    label 'process_medium'
    label 'process_long'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', saveAs: { filename -> "price/${meta.id}/${filename}" }
    container 'ghcr.io/jackcurragh/translon-gedi-price:1.0.0'

    input:
    tuple val(meta), val(contig), path(bam), path(bai), path(gtf), path(fasta), path(index)

    output:
    tuple val(meta), path("${prefix}.orfs.tsv"), emit: orfs_tsv
    tuple val(meta), path("${prefix}.orfs.cit"), optional: true, emit: orfs_cit
    tuple val(meta), path("${prefix}.orfs.cit.metadata.json"), optional: true, emit: orfs_metadata
    tuple val(meta), path("${prefix}.codons.cit"), optional: true, emit: codons_cit
    tuple val(meta), path("${prefix}.model"), optional: true, emit: model
    tuple val(meta), path("${prefix}.signal.tsv"), optional: true, emit: signal
    tuple val(meta), path("${prefix}.param"), optional: true, emit: param

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${contig}"
    def threads = task.cpus ?: 1
    """
    mkdir -p prepared_bams
    reference_fasta=\$(ls ${index}/*.fa ${index}/*.fasta ${index}/*.fna 2>/dev/null | head -n 1)
    test -n "\$reference_fasta" || { echo 'GEDI PRICE index does not contain its reference FASTA' >&2; exit 1; }
    samtools calmd -b -@ ${threads} ${bam} \"\$reference_fasta\" > prepared_bams/input.bam
    samtools index -@ ${threads} prepared_bams/input.bam
    gedi -e Bam2CIT -p price_input.bamlist.cit prepared_bams/input.bam
    sed "s|file=\\\"[^\\\"]*/|file=\\\"\$PWD/${index}/|g" ${index}/${contig}.oml > ${prefix}.genomic.oml
    gedi -e Price \\
        -reads price_input.bamlist.cit \\
        -genomic ${prefix}.genomic.oml \\
        -prefix ${prefix} \\
        -nthreads ${threads}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${contig}"
    """
    touch ${prefix}.orfs.tsv
    """
}
