process CHOROS {
    tag "${meta.id}"
    label 'process_high'

    container "${params.choros_container}"

    publishDir "${params.outdir}/choros", mode: 'copy', pattern: "*.choros_*"

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bai), path(best_offsets)
    path ribometric_annotation
    path transcriptome_fasta

    output:
    tuple val(meta), path("*.choros_counts.tsv.gz"), emit: corrected_counts
    tuple val(meta), path("*.choros_coefficients.tsv"), emit: coefficients
    tuple val(meta), path("*.choros_metrics.tsv"), emit: metrics
    tuple val(meta), path("*.choros_bam_metrics.tsv"), emit: bam_metrics
    tuple val(meta), path("*.choros_offsets.tsv"), emit: offset_rules
    path "versions.yml", emit: versions

    when:
    params.run_choros

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes = params.choros_num_genes != null ? params.choros_num_genes : 250
    def minCoverage = params.choros_min_coverage != null ? params.choros_min_coverage : 5
    def minNonzero = params.choros_min_nonzero != null ? params.choros_min_nonzero : 100
    """
    samtools collate -@ $task.cpus -o ${prefix}.collated.bam ${transcriptome_bam}

    python3 $projectDir/bin/prepare_choros_bam.py \
        --bam ${prefix}.collated.bam \
        --output ${prefix}.choros_input.bam \
        --metrics ${prefix}.choros_bam_metrics.tsv

    python3 $projectDir/bin/prepare_choros_inputs.py \
        --annotation ${ribometric_annotation} \
        --offsets ${best_offsets} \
        --lengths-output ${prefix}.choros_lengths.tsv \
        --offsets-output ${prefix}.choros_offsets.tsv

    Rscript $projectDir/bin/run_choros.R \
        ${prefix}.choros_input.bam \
        ${transcriptome_fasta} \
        ${prefix}.choros_lengths.tsv \
        ${prefix}.choros_offsets.tsv \
        ${prefix} \
        ${genes} \
        ${minCoverage} \
        ${minNonzero} \
        $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        choros: \$(Rscript -e 'cat(as.character(packageVersion("choros")))')
        samtools: \$(samtools --version | sed -n '1s/samtools //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e 'transcript\tcod_idx\tcount\tcorrected' | gzip -c > ${prefix}.choros_counts.tsv.gz
    echo -e 'group\tterm\testimate' > ${prefix}.choros_coefficients.tsv
    echo -e 'metric\tvalue\ntraining_transcripts\t250' > ${prefix}.choros_metrics.tsv
    echo -e 'metric\tvalue\nretained_fraction\t1.0' > ${prefix}.choros_bam_metrics.tsv
    echo -e 'length\tframe_0\tframe_1\tframe_2\n28\t15\t14\t16' > ${prefix}.choros_offsets.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        choros: pilot
        samtools: 1.20
    END_VERSIONS
    """
}
