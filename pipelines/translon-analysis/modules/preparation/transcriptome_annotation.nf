process MAKE_TRANSCRIPTOME_ANNOTATION {
    tag "${meta.id}"
    label 'process_medium'
    container 'python:3.12-bookworm'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path(bam), path(bai), path('transcriptome.gtf'), path('transcriptome.fa'), emit: annotation
    script:
    """
    python3 ${params.translon_analysis_bin}/make_transcriptome_annotation.py --gtf ${gtf} --fasta ${fasta} --out-gtf transcriptome.gtf --out-fasta transcriptome.fa
    """
    stub:
    """
    touch transcriptome.gtf transcriptome.fa
    """
}
