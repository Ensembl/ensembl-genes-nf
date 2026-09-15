process EXTRACT_COMBINED_TRANSCRIPTS {
    tag "${meta.id}:transcripts"
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1'

    input:
    tuple val(meta), path(combined_bed)
    path reference

    output:
    tuple val(meta), path('combined_transcripts.fa'), emit: transcripts
    path 'sequence_provenance.tsv', emit: provenance
    path 'versions.yml', emit: versions

    script:
    """
    bedtools getfasta -fi ${reference} -bed ${combined_bed} -name -split -s -fo combined_transcripts.fa
    printf 'reference_fasta\t%s\nmodels_extracted\t%s\n' \\
        "\$(sha256sum ${reference} | awk '{print \$1}')" "\$(grep -c '^>' combined_transcripts.fa || true)" > sequence_provenance.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    """
    printf '>stub.1\nATGAAATAA\n' > combined_transcripts.fa
    printf 'reference_fasta\tstub\nmodels_extracted\t1\n' > sequence_provenance.tsv
    printf '"%s":\n    bedtools: 2.31.1\n' '${task.process}' > versions.yml
    """
}
