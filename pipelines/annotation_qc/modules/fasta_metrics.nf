process FASTA_METRICS {
    tag "$gca_accession"

    publishDir params.outdir, mode: 'copy', saveAs: { filename ->
        def prefix = "${gca_accession}."
        "${gca_accession}/${filename.startsWith(prefix) ? filename.substring(prefix.length()) : filename}"
    }

    input:
	tuple val(gca_accession), path(normalised), path(fasta)

	output:
	tuple val(gca_accession), path("${gca_accession}.fasta_metrics.csv"), path("${gca_accession}.masked_overlaps.tsv.gz"), emit: fasta_metrics

	script:
	"""
    python ${projectDir}/bin/fasta_metrics.py ${normalised} ${fasta} \\
        --output-csv ${gca_accession}.fasta_metrics.csv \\
        --masked-overlaps-output ${gca_accession}.masked_overlaps.tsv.gz
	"""
}
