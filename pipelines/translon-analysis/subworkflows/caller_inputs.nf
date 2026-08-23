/* Convert the published RiboSeq alignment contract only where a caller needs
 * a legacy representation. Direct BAM consumers are passed through unchanged. */

process PREPARE_TRANSCRIPT_MODELS {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' }
    // Nextflow's task-metrics wrapper requires `ps`; the slim image omits
    // procps, so use the full Debian runtime for this process.
    container 'python:3.12-bookworm'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    output:
    tuple val(meta), path(bam), path(bai), path('models.genePred'), path('models.bed12'), emit: models
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    python3 ${projectDir}/bin/gtf_to_transcript_models.py \
        --gtf ${gtf} --genepred models.genePred --bed12 models.bed12
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcript_models: 1
    END_VERSIONS
    """
}

process PREPARE_RIBORF_READS {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' }
    // This step uses both samtools and RibORF's offsetCorrect.pl. The
    // RibORF image provides samtools plus RIBORF_HOME; the samtools-only
    // image cannot satisfy the offset-correction contract.
    container 'ghcr.io/jackcurragh/translon-riborf:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai), path(genepred), path(bed12), path(offsets)
    output:
    tuple val(meta), path('reads.sam'), path(genepred), path(bed12), emit: riborf
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    samtools view -h ${bam} > uncorrected.sam
    perl \$RIBORF_HOME/offsetCorrect.pl -r uncorrected.sam -p ${offsets} -o reads.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch reads.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}

workflow PREPARE_CALLER_INPUTS {
    take:
    transcriptome
    gtf
    offsets

    main:
    PREPARE_TRANSCRIPT_MODELS(transcriptome, gtf)
    models_with_offsets = PREPARE_TRANSCRIPT_MODELS.out.models
        .map { meta, bam, bai, genepred, bed12 -> tuple(meta.id, meta, bam, bai, genepred, bed12) }
        .join(offsets.map { meta, offset -> tuple(meta.id, offset) }, by: 0)
        .map { id, meta, bam, bai, genepred, bed12, offset -> tuple(meta, bam, bai, genepred, bed12, offset) }
    PREPARE_RIBORF_READS(models_with_offsets)

    emit:
    transcript_models = PREPARE_TRANSCRIPT_MODELS.out.models
    riborf = PREPARE_RIBORF_READS.out.riborf
}
