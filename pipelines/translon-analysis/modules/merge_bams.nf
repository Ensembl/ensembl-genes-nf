/* Merge only like-for-like alignments.  Transcriptome and genome BAMs are
 * deliberately merged on separate channels so coordinate systems can never be
 * mixed accidentally. */

process MERGE_RIBO_BAMS {
    tag "${meta.merge_group}:${meta.bam_type}"
    label 'process_medium'
    errorStrategy 'terminate'
    container params.container_samtools
    publishDir "${params.outdir}/merged_inputs", mode: 'copy', saveAs: { filename -> "${meta.merge_group}/${meta.bam_type}/${filename}" }

    input:
    tuple val(meta), path(bams, stageAs: 'bams/*'), path(bais, stageAs: 'bais/*')

    output:
    tuple val(meta), path('merged.bam'), path('merged.bam.bai'), emit: merged
    path 'merge_manifest.tsv', emit: manifest

    script:
    """
    set -euo pipefail
    bam_count=\$(find bams -maxdepth 1 -type f -name '*.bam' | wc -l | tr -d ' ')
    test "\$bam_count" -gt 0 || { echo 'No BAMs supplied to MERGE_RIBO_BAMS' >&2; exit 1; }
    samtools merge -f -@ ${task.cpus ?: 4} merged.bam bams/*.bam
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

process MERGE_RIBO_OFFSETS {
    tag "${meta.merge_group}:offsets"
    label 'process_low'
    publishDir "${params.outdir}/merged_inputs", mode: 'copy', saveAs: { filename -> "${meta.merge_group}/offsets/${filename}" }

    input:
    tuple val(meta), path(offsets, stageAs: 'offsets/*')

    output:
    tuple val(meta), path('merged.offsets.tsv'), emit: offsets

    script:
    """
    set -euo pipefail
    first=\$(find offsets -maxdepth 1 -type f | sort | head -n1)
    test -n "\$first" || { echo 'No offset files supplied' >&2; exit 1; }
    head -n1 "\$first" > merged.offsets.tsv
    : > merged.offsets.body.tsv
    for offset_file in offsets/*; do
        tail -n +2 "\$offset_file" | awk -F '\\t' '{ print \$1 "\\t" \$2 }' >> merged.offsets.body.tsv
    done
    awk -F '\\t' '
        NF >= 2 {
            if ((\$1 in seen) && seen[\$1] != \$2) {
                printf "Conflicting offsets for read length %s: %s vs %s\\n", \$1, seen[\$1], \$2 > "/dev/stderr"
                bad = 1
            }
            seen[\$1] = \$2
        }
        END {
            if (bad) exit 2
            for (length in seen) print length "\\t" seen[length]
        }
    ' merged.offsets.body.tsv > merged.offsets.unique.tsv
    sort -k1,1n merged.offsets.unique.tsv >> merged.offsets.tsv
    rm merged.offsets.body.tsv merged.offsets.unique.tsv
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
        offset_groups = offsets
            .map { meta, offset -> tuple([id: meta.merge_group ?: meta.id, merge_group: meta.merge_group ?: meta.id], offset) }
            .groupTuple(by: 0)
            .map { meta, files -> tuple(meta, files) }
        MERGE_RIBO_OFFSETS(offset_groups)
        offsets_out = MERGE_RIBO_OFFSETS.out.offsets
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
