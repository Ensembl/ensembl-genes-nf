/* Merge only like-for-like alignments.  Transcriptome and genome BAMs are
 * deliberately merged on separate channels so coordinate systems can never be
 * mixed accidentally. */

process MERGE_RIBO_BAMS {
    tag "${meta.merge_group}:${meta.bam_type}"
    label 'process_medium'
    errorStrategy 'terminate'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "${params.outdir}/merged_inputs", mode: 'copy', saveAs: { filename -> "${meta.merge_group}/${meta.bam_type}/${filename}" }

    input:
    tuple val(meta), path(bams, stageAs: 'bams/*'), path(bais, stageAs: 'bais/*')

    output:
    tuple val(meta), path('merged.bam'), path('merged.bam.bai'), emit: merged
    path 'merge_manifest.tsv', emit: manifest

    script:
    """
    set -euo pipefail
    bam_files=(bams/*.bam)
    test -e "\${bam_files[0]}" || { echo 'No BAMs supplied to MERGE_RIBO_BAMS' >&2; exit 1; }
    samtools merge -f -@ ${task.cpus ?: 4} merged.bam "\${bam_files[@]}"
    samtools index -@ ${task.cpus ?: 4} merged.bam
    {
        printf 'merge_group\\tbam_type\\tmerged_bam\\tsource_bam\\n'
        for bam in bams/*.bam; do
            printf '%s\\t%s\\t%s\\t%s\\n' '${meta.merge_group}' '${meta.bam_type}' 'merged.bam' "\$bam"
        done
    } > merge_manifest.tsv
    """

    stub:
    """
    touch merged.bam merged.bam.bai
    printf 'merge_group\\tbam_type\\tmerged_bam\\tsource_bam\\n' > merge_manifest.tsv
    """
}

workflow MERGE_RIBO_INPUTS {
    take:
    transcriptome
    genome
    offsets
    merge_inputs
    merge_group

    main:
    if (merge_inputs) {
        tx_grouped = transcriptome
            .map { meta, bam, bai -> tuple([id: meta.merge_group ?: merge_group, merge_group: meta.merge_group ?: merge_group, bam_type: 'transcriptome'], bam, bai) }
            .groupTuple(by: 0)
            .map { meta, bams, bais -> tuple(meta, bams, bais) }
        gn_grouped = genome
            .map { meta, bam, bai -> tuple([id: meta.merge_group ?: merge_group, merge_group: meta.merge_group ?: merge_group, bam_type: 'genome'], bam, bai) }
            .groupTuple(by: 0)
            .map { meta, bams, bais -> tuple(meta, bams, bais) }
        MERGE_RIBO_BAMS(tx_grouped.mix(gn_grouped))
        tx_out = MERGE_RIBO_BAMS.out.merged.filter { meta, bam, bai -> meta.bam_type == 'transcriptome' }
        gn_out = MERGE_RIBO_BAMS.out.merged.filter { meta, bam, bai -> meta.bam_type == 'genome' }
        // Per-sample offsets must not be pooled.  They describe the individual
        // library's read-length/P-site behaviour.  The parent workflow runs a
        // fresh RiboMetric calculation on the pooled transcriptome BAM instead.
        offsets_out = channel.empty()
    } else {
        tx_out = transcriptome
        gn_out = genome
        offsets_out = offsets
    }

    emit:
    transcriptome = tx_out
    genome = gn_out
    offsets = offsets_out
    manifests = merge_inputs ? MERGE_RIBO_BAMS.out.manifest : channel.empty()
}
