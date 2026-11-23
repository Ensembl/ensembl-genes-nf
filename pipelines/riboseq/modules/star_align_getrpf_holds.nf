/* Diagnostic STAR alignment for raw reads whose getRPF transform was withheld.
 *
 * This process is deliberately separate from the production alignment. It
 * always uses Local ends so the resulting CIGAR strings can be inspected for
 * protocol-specific soft-clipping evidence. Its BAMs must not be interpreted
 * as production RPF alignments or sent to RiboMetric/QC gating unchanged.
 */

process STAR_ALIGN_GETRPF_HOLDS {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::star=2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_1' :
        'biocontainers/star:2.7.11b--h43eeafb_1' }"

    input:
    tuple val(meta), path(reads)
    path index
    path gtf

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam", arity: '1'), emit: bam
    tuple val(meta), path("*.Aligned.toTranscriptome.out.bam", arity: '1'), emit: transcriptome_bam
    tuple val(meta), path("*.Log.final.out", arity: '1'), emit: log
    tuple val(meta), path("*.Log.out", arity: '1'), emit: log_out
    tuple val(meta), path("*.Log.progress.out", arity: '1'), emit: log_progress
    tuple val(meta), path("*.SJ.out.tab", arity: '1'), emit: splice_junctions
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def allow_introns = params.allow_introns ? '--alignIntronMax 1000000 --alignMatesGapMax 1000000' : ''
    def unzip_command = reads.name.endsWith('.gz') ? 'zcat' : 'cat'
    def soft_clip_at_ref = params.align_soft_clip_at_ref_ends ? "--alignSoftClipAtReferenceEnds ${params.align_soft_clip_at_ref_ends}" : ''

    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM \\
        --outFilterMultimapNmax ${params.max_multimappers} \\
        --outFilterMismatchNmax ${params.mismatches} \\
        --readFilesCommand $unzip_command \\
        --quantMode TranscriptomeSAM \\
        --alignEndsType Local \\
        $allow_introns \\
        $soft_clip_at_ref \\
        $args \\
        --outFileNamePrefix ${prefix}. \\
        --sjdbGTFfile $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Aligned.toTranscriptome.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.SJ.out.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11b
    END_VERSIONS
    """
}
