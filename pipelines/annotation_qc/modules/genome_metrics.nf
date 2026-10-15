process GENOME_METRICS {
    tag "$gca_accession"

    publishDir params.outdir, mode: 'copy', saveAs: { filename ->
        def prefix = "${gca_accession}."
        "${gca_accession}/${filename.startsWith(prefix) ? filename.substring(prefix.length()) : filename}"
    }

    input:
	tuple val(gca_accession), path(normalised), path(fasta_metrics)

	output:
	tuple val(gca_accession), path("${gca_accession}.genome_metrics.csv"), path("${gca_accession}.feature_overlaps.tsv.gz"), emit: genome_metrics

	script:
	"""
    python ${projectDir}/bin/genome_metrics.py ${normalised} --feature-type gene \\
        --fasta-metrics ${fasta_metrics} \\
        --output-csv ${gca_accession}.genome_metrics.csv \\
        --feature-overlaps-output ${gca_accession}.feature_overlaps.tsv.gz
	"""
}
