/*
 * INPUT CONTRACT
 *
 * Resolve the published output tree of pipelines/riboseq into the channel
 * contracts expected by the ORF callers. Explicit globs override the
 * conventional output-tree discovery rules.
 */

workflow LOAD_RIBOSEQ_OUTPUTS {
    take:
    riboseq_outdir
    transcriptome_bam_glob
    genome_bam_glob
    offsets_glob
    translonscorer_glob

    main:
    if (!riboseq_outdir) {
        error 'riboseq_outdir is required: point this pipeline at pipelines/riboseq --outdir'
    }

    def root = file(riboseq_outdir, checkIfExists: true)
    def tx_pattern = transcriptome_bam_glob ?: "${root}/**/*Aligned.toTranscriptome.out.bam"
    def gn_pattern = genome_bam_glob ?: "${root}/**/*Aligned.sortedByCoord.out.bam"

    transcriptome = channel.fromPath(tx_pattern, checkIfExists: true)
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai")
            def id = bam.baseName.replaceAll(/\.Aligned\.toTranscriptome\.out$/, '')
            tuple([id: id, bam_type: 'transcriptome'], bam, bai)
        }

    genome = channel.fromPath(gn_pattern, checkIfExists: true)
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai")
            def id = bam.baseName.replaceAll(/\.Aligned\.sortedByCoord\.out$/, '')
            tuple([id: id, bam_type: 'genome'], bam, bai)
        }

    offsets = channel.fromPath(
        offsets_glob ?: "${root}/**/*.offsets.{pass,selected,good,great}.tsv",
        checkIfExists: false
    )
    translonscorer = channel.fromPath(
        translonscorer_glob ?: "${root}/**/*_orfs_scored.csv",
        checkIfExists: false
    )

    emit:
    transcriptome = transcriptome
    genome = genome
    offsets = offsets
    translonscorer = translonscorer
}
