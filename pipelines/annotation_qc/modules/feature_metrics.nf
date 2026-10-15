process FEATURE_METRICS {
    tag "$gca_accession"

    publishDir params.outdir, mode: 'copy', saveAs: { filename ->
        def prefix = "${gca_accession}."
        "${gca_accession}/${filename.startsWith(prefix) ? filename.substring(prefix.length()) : filename}"
    }

    input:
	tuple val(gca_accession), path(normalised)

	output:
	tuple val(gca_accession), path("${gca_accession}.feature_metrics.csv"), path("${gca_accession}.flagged_features.tsv.gz"), emit: feature_metrics

	script:
	"""
	python ${projectDir}/bin/feature_metrics.py ${normalised} \\
        --output-csv ${gca_accession}.feature_metrics.csv \\
        --flagged-features ${gca_accession}.flagged_features.tsv.gz
	"""
}
