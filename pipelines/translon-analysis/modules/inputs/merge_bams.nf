process MERGE_RIBO_BAMS {
    tag "${meta.merge_group}:${meta.bam_type}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    input:
    tuple val(meta), path(bams, stageAs: 'bams/*'), path(bais, stageAs: 'bais/*')
    output:
    tuple val(meta), path('merged.bam'), path('merged.bam.bai'), emit: merged
    path 'merge_manifest.tsv', emit: manifest
    script:
    """
    set -euo pipefail
    bam_files=(bams/*.bam)
    test -e "\${bam_files[0]}" || { echo 'No BAMs supplied' >&2; exit 1; }
    samtools merge -f -@ ${task.cpus ?: 4} merged.bam "\${bam_files[@]}"
    samtools index -@ ${task.cpus ?: 4} merged.bam
    printf 'merge_group\tbam_type\tmerged_bam\tsource_bam\n' > merge_manifest.tsv
    """
    stub:
    """
    touch merged.bam merged.bam.bai
    printf 'merge_group\tbam_type\tmerged_bam\tsource_bam\n' > merge_manifest.tsv
    """
}
