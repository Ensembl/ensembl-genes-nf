process MAKE_PARTITION_MANIFEST {
    tag "${meta.id}:${params.partition_mode ?: 'transcriptome'}"
    label 'process_low'
    // The manifest estimates BAM load with samtools; use the same pinned
    // samtools image used by BAM preparation rather than a Python-only image.
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    input:
    tuple val(meta), path(annotation), path(bam)
    output:
    tuple val(meta), path('partition_manifest.tsv'), emit: manifest
    script:
    def mode = params.partition_mode ?: 'transcriptome'
    def bamArg = "--bam ${bam}"
    def annotationArg = mode == 'transcriptome' ? "--gtf ${annotation}" : ''
    def faiArg = mode == 'genome' ? "--fai ${annotation}" : ''
    """
    python3 ${params.translon_analysis_bin}/make_partition_manifest.py --mode ${mode} ${annotationArg} ${bamArg} ${faiArg} --partitions ${params.partition_count ?: 1} --window-size ${params.partition_window_size ?: 0} --padding ${params.partition_padding ?: 0} --output partition_manifest.tsv
    """
    stub:
    """
    printf 'partition_id\tmode\tcontig\tstart\tend\tpadding\tannotation_load\testimated_read_load\n0\t${params.partition_mode ?: 'transcriptome'}\tchr1\t0\t100\t0\t1\t\n' > partition_manifest.tsv
    """
}
