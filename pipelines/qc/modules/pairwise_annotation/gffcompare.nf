process PAIRWISE_GFFCOMPARE {
    tag { "${meta.id}_${direction}" }
    label 'process_low'

    container "quay.io/biocontainers/gffcompare:0.12.6--h4ac6f70_2"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*"

    input:
        tuple val(meta), val(direction), path(ref_gtf, stageAs: 'ref.gtf'), path(qry_gtf, stageAs: 'qry.gtf')

    output:
        tuple val(meta), val(direction), path("tmap_${direction}.tsv"), emit: tmap
        tuple val(meta), val(direction), path("stats_${direction}.txt"), emit: stats
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail

        gffcompare -r ref.gtf -o ${direction} qry.gtf || true
        cp ${direction}.qry.gtf.tmap tmap_${direction}.tsv 2>/dev/null || : > tmap_${direction}.tsv
        cp ${direction}.stats stats_${direction}.txt 2>/dev/null || : > stats_${direction}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffcompare: \$(gffcompare --version 2>&1 | sed 's/gffcompare v//g')
        END_VERSIONS
        """

    stub:
        """
        touch tmap_${direction}.tsv
        touch stats_${direction}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffcompare: 0.12.6
        END_VERSIONS
        """
}

process PAIRWISE_PARSE_GFFCOMPARE {
    tag { "${meta.id}_${direction}" }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "class_counts_*"

    input:
        tuple val(meta), val(direction), path(tmap)

    output:
        tuple val(meta), val(direction), path("class_counts_${direction}.tsv"), emit: counts
        path "versions.yml", emit: versions

    script:
        """
        set -euo pipefail

        parse_gffcompare_tmap.py \\
            --tmap ${tmap} \\
            --assembly-accession ${meta.id} \\
            --sample-name ${meta.id} \\
            --direction ${direction} \\
            --output class_counts_${direction}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tdirection\\tclass_code\\tn_transcripts\\tdenominator\\tpct\\n" > class_counts_${direction}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
