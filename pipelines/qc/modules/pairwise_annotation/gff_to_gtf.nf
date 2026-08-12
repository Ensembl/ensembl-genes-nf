process PAIRWISE_GFF_TO_GTF {
    tag { "${meta.id}_${source_label}" }
    label 'process_medium'

    container "quay.io/biocontainers/gffread:0.12.7--hd03093a_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*"

    input:
        tuple val(meta), val(source_label), path(gff)

    output:
        tuple val(meta), val(source_label), path("${meta.id}.${source_label}.gtf"), emit: gtf
        path "versions.yml", emit: versions

    script:
        """
        set -euo pipefail

        gff_to_gtf.sh ${gff} ${meta.id}.${source_label}.gtf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffread: \$(gffread --version 2>&1 | sed 's/^gffread //g')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.${source_label}.gtf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffread: 0.12.7
        END_VERSIONS
        """
}
