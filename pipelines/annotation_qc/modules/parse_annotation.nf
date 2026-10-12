process PARSE_ANNOTATION {
    tag "$gca_accession"

	input:
	tuple val(gca_accession), path(annotation)
	
	output:
	tuple val(gca_accession), path("${gca_accesion}.normalised.tsv.gz"), emit:
	normalised

	script:
	"""
	python parse_annotation.py ${annotation} \\
		--output ${gca_accession}.normalised.tsv.gz
	"""

}
