// MERGE_RNASEQ
// Merge per-sample RNA-seq GFF3 files into a single non-redundant set.
// Clusters overlapping transcripts and assigns unified gene IDs.

process MERGE_RNASEQ {
    label 'process_medium'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/rnaseq", mode: 'copy', pattern: "merged_rnaseq.gff3"

    input:
    path gff3_files   // list of per-sample GFF3 files

    output:
    path "merged_rnaseq.gff3", emit: gff3
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    merge_rnaseq_gff3.py \\
        --gff3 ${gff3_files} \\
        --out  merged_rnaseq.gff3 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > merged_rnaseq.gff3
    printf 'chr1\\tStringTie2\\tgene\\t1000\\t5000\\t.\\t+\\t.\\tID=rnaseq_merged_gene_00000001;biotype=rnaseq_merged\\n' \\
        >> merged_rnaseq.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
