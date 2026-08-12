process ANNOTATION_INTEGRATE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy'

    // Wave resolves this environment to a reproducible container on HPC.
    // The Python application itself uses only the standard library.
    conda 'conda-forge::python=3.12 bioconda::gffcompare=0.12.6 bioconda::genometools-genometools=1.6.5'

    input:
    tuple val(meta), path(projected_gff), path(manual_gff), path(decision_tsv), path(assembly_fasta), path(assembly_fai)
    path integration_script

    output:
    tuple val(meta), path("${meta.id}.annotation_review.tsv"), emit: review
    tuple val(meta), path("${meta.id}.review_cases.csv"), emit: cases
    tuple val(meta), path("${meta.id}.review_cases_representative.csv"), emit: representative_cases
    tuple val(meta), path("${meta.id}.havana_decisions.tsv"), emit: havana_decisions
    tuple val(meta), path("${meta.id}.integrated.gff3"), emit: integrated
    tuple val(meta), path("${meta.id}.jbrowse_session.json"), emit: session
    path "${meta.id}.jbrowse_config.json", emit: jbrowse_config
    path "${meta.id}.gff_validation.*.log", emit: validation_logs
    path "${meta.id}.manual_vs_projected.*", optional: true, emit: gffcompare
    path "${meta.id}.README.txt", emit: readme

    script:
    prefix = meta.id
    def projected_arg = projected_gff.name.startsWith('empty_') ? '' : "--projected-gff ${projected_gff}"
    def manual_arg = manual_gff.name.startsWith('empty_') ? '' : "--manual-gff ${manual_gff}"
    def decision_arg = decision_tsv.name.startsWith('empty_') ? '' : "--decision-tsv ${decision_tsv}"
    """
    set -euo pipefail

    mkdir -p integration_out

    if ${params.run_external_validation} && [ -s ${projected_gff} ] && grep -qv '^#' ${projected_gff}; then
        gt gff3validator -strict ${projected_gff} > ${prefix}.gff_validation.projected.log 2>&1
    else
        echo 'projected annotation absent' > ${prefix}.gff_validation.projected.log
    fi

    if ${params.run_external_validation} && [ -s ${manual_gff} ] && grep -qv '^#' ${manual_gff}; then
        gt gff3validator -strict ${manual_gff} > ${prefix}.gff_validation.manual.log 2>&1
    else
        echo 'manual annotation absent' > ${prefix}.gff_validation.manual.log
    fi

    if ${params.run_external_validation} && [ -s ${projected_gff} ] && grep -qv '^#' ${projected_gff} && [ -s ${manual_gff} ] && grep -qv '^#' ${manual_gff}; then
        gffcompare -r ${projected_gff} -o ${prefix}.manual_vs_projected ${manual_gff}
    fi

    python3 ${integration_script} \\
        --assembly ${params.assembly_name} \\
        --output-dir integration_out \\
        ${projected_arg} ${manual_arg} ${decision_arg} \\
        --fasta ${assembly_fasta} --fai ${assembly_fai}

    mv integration_out/annotation_review.tsv ${prefix}.annotation_review.tsv
    mv integration_out/review_cases.csv ${prefix}.review_cases.csv
    mv integration_out/review_cases_representative.csv ${prefix}.review_cases_representative.csv
    mv integration_out/havana_decisions.tsv ${prefix}.havana_decisions.tsv
    mv integration_out/integrated.gff3 ${prefix}.integrated.gff3
    mv integration_out/jbrowse_session.json ${prefix}.jbrowse_session.json
    mv integration_out/jbrowse_config.json ${prefix}.jbrowse_config.json
    mv integration_out/README.txt ${prefix}.README.txt
    """

    stub:
    """
    touch ${meta.id}.annotation_review.tsv
    touch ${meta.id}.review_cases.csv
    touch ${meta.id}.review_cases_representative.csv
    touch ${meta.id}.havana_decisions.tsv
    touch ${meta.id}.integrated.gff3
    touch ${meta.id}.jbrowse_session.json
    touch ${meta.id}.jbrowse_config.json
    touch ${meta.id}.gff_validation.projected.log
    touch ${meta.id}.gff_validation.manual.log
    touch ${meta.id}.README.txt
    """
}
