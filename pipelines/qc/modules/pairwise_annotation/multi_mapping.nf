process PAIRWISE_MULTI_MAPPING {
    tag { meta.id }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*_multi_mapping.tsv"

    input:
        tuple val(meta), path(all_pairs)

    output:
        tuple val(meta), path("${meta.id}_multi_mapping.tsv"), emit: metrics
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def minOverlap = params.pairwise_multi_mapping_min_overlap ?: 0.1
        """
        set -euo pipefail

        analyze_multi_mapping.py \\
            --all-pairs ${all_pairs} \\
            --output ${meta.id}_multi_mapping.tsv \\
            --assembly-accession ${meta.id} \\
            --sample-name ${meta.id} \\
            --min-overlap ${minOverlap}

        normalise_pairwise_tsv.py --input ${meta.id}_multi_mapping.tsv --output ${meta.id}_multi_mapping.normalised.tsv
        mv ${meta.id}_multi_mapping.normalised.tsv ${meta.id}_multi_mapping.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tgene_id\\tsource\\tn_matches\\trelationship\\tclassification\\n" > ${meta.id}_multi_mapping.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
