process REPORT_CONSENSUS {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pip_duckdb_pandas:57c3741f55d53490"
    
    tag "${meta.id}"

    publishDir "${params.outdir}/consensus_reports/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(samplesheet), path(bed_files)
    val ucsc_session_url
    path gencode_gtf, stageAs: 'gencode.gtf*'

    output:
    tuple val(meta), path("*.tsv"),           emit: results

    script:
    def gtf_arg = gencode_gtf.name != 'NO_FILE' ? "-g ${gencode_gtf}" : ""
    """
    consensus.py \\
        -s ${samplesheet} \\
        -n ${meta.id} \\
        -o . \\
        -u "${ucsc_session_url}" \\
        ${gtf_arg}
    """
    
    stub:
    """
    mkdir -p consensus_results
    touch consensus_results/${prefix}.txt
    """
}
