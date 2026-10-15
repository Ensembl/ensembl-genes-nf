process PARSE_ANNOTATION {
    tag "$gca_accession"

    publishDir params.outdir, mode: 'copy', saveAs: { filename ->
        def prefix = "${gca_accession}."
        "${gca_accession}/${filename.startsWith(prefix) ? filename.substring(prefix.length()) : filename}"
    }

    input:
	tuple val(gca_accession), file(annotation), file(fasta)

	output:
	tuple val(gca_accession), path("${gca_accession}.normalised.tsv.gz"), path(fasta), emit: parsed

	script:
	"""
	python ${projectDir}/bin/parse_annotation.py ${annotation} \\
		--output ${gca_accession}.normalised.tsv.gz
	"""

}
