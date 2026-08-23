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
    samplesheet
    transcriptome_bam_glob
    genome_bam_glob
    offsets_glob
    translonscorer_glob
    merge_group

    main:
    if (!riboseq_outdir && !samplesheet) {
        error 'Provide either --riboseq_outdir or --samplesheet'
    }

    def root = riboseq_outdir ? file(riboseq_outdir, checkIfExists: true) : null
    def discovered_bams = []
    if (root) root.toFile().eachFileRecurse { candidate ->
        if (candidate.isFile() && candidate.name.endsWith('.bam')) discovered_bams << candidate.toPath()
    }
    def tx_pattern = transcriptome_bam_glob
    def gn_pattern = genome_bam_glob
    def sheet_rows = samplesheet ? channel.fromPath(samplesheet, checkIfExists: true).splitCsv(header: true, sep: '\t') : null

    transcriptome = sheet_rows
        ? sheet_rows.filter { row -> row.transcriptome_bam }
            .map { row ->
                def bam = file(row.transcriptome_bam, checkIfExists: true)
                def bai = row.transcriptome_bai ? file(row.transcriptome_bai, checkIfExists: true) : file("${bam}.bai", checkIfExists: true)
                tuple([id: row.sample_id, merge_group: row.merge_group ?: row.sample_id, bam_type: 'transcriptome'], bam, bai)
            }
        : (tx_pattern
            ? channel.fromPath(tx_pattern, checkIfExists: true)
            : channel.fromList(discovered_bams)
                .filter { bam -> bam.name.endsWith('Aligned.toTranscriptome.out.bam') })
        .map { bam ->
            def bai = file("${bam}.bai", checkIfExists: false)
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai", checkIfExists: false)
            def id = bam.baseName.replaceAll(/\.Aligned\.toTranscriptome\.out$/, '')
            tuple([id: id, merge_group: merge_group, bam_type: 'transcriptome'], bam, bai)
        }

    genome = sheet_rows
        ? sheet_rows.filter { row -> row.genome_bam }
            .map { row ->
                def bam = file(row.genome_bam, checkIfExists: true)
                def bai = row.genome_bai ? file(row.genome_bai, checkIfExists: true) : file("${bam}.bai", checkIfExists: true)
                tuple([id: row.sample_id, merge_group: row.merge_group ?: row.sample_id, bam_type: 'genome'], bam, bai)
            }
        : (gn_pattern
            ? channel.fromPath(gn_pattern, checkIfExists: true)
            : channel.fromList(discovered_bams)
                .filter { bam -> bam.name.endsWith('Aligned.sortedByCoord.out.bam') })
        .map { bam ->
            def bai = file("${bam}.bai", checkIfExists: false)
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai", checkIfExists: false)
            def id = bam.baseName.replaceAll(/\.Aligned\.sortedByCoord\.out$/, '')
            tuple([id: id, merge_group: merge_group, bam_type: 'genome'], bam, bai)
        }

    ribo_fastq = sheet_rows
        ? sheet_rows.filter { row -> row.ribo_fastq }
            .map { row -> tuple([id: row.sample_id, bam_type: 'ribo_fastq'], file(row.ribo_fastq, checkIfExists: true)) }
        : channel.empty()

    offsets = offsets_glob
        ? channel.fromPath(offsets_glob, checkIfExists: false)
        : (sheet_rows
            ? sheet_rows.filter { row -> row.offsets }
                .map { row -> tuple([id: row.sample_id, merge_group: row.merge_group ?: row.sample_id], file(row.offsets, checkIfExists: true)) }
            : (root ? channel.fromPath("${root}/**/*.offsets.{pass,selected,good,great}.tsv", checkIfExists: false) : channel.empty()))
    translonscorer = translonscorer_glob
        ? channel.fromPath(translonscorer_glob, checkIfExists: false)
        : (root ? channel.fromPath("${root}/**/*_orfs_scored.csv", checkIfExists: false) : channel.empty())

    emit:
    transcriptome = transcriptome
    genome = genome
    ribo_fastq = ribo_fastq
    offsets = offsets
    translonscorer = translonscorer
}
