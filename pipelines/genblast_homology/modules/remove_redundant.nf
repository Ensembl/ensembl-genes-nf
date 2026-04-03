// REMOVE_REDUNDANT
// Merge GFF3 files from all GenBlast batches, cluster overlapping gene models,
// and apply biotype-priority filtering: within each overlap cluster retain only
// models of the highest priority tier (genblast_1 > ... > genblast_7).
// Emulates Ensembl RemoveRedundantGenes layer-annotation approach.

process REMOVE_REDUNDANT {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/genblast_homology", mode: 'copy', pattern: "*.final.gff3"

    input:
    tuple val(meta), path(gff3_files)   // all classified GFF3s collected

    output:
    tuple val(meta), path("*.final.gff3"), emit: gff3
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    echo '##gff-version 3' > merged.gff3
    cat ${gff3_files} | grep -v '^#' >> merged.gff3

    remove_redundant.py \\
        --gff3      merged.gff3 \\
        --out       ${prefix}.final.gff3 \\
        --sample-id ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\nchr1\tgenBlastG\tgene\t1000\t2000\t.\t+\t.\tID=${prefix}_gbh_gene_00000001;Name=STUB;biotype=genblast_1\n' \\
        > ${prefix}.final.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
