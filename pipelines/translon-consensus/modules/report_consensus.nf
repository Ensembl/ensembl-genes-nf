// TODO: Rename this process to match your tool/module name (e.g., MY_TOOL)
process REPORT_CONSENSUS {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pip_duckdb_pandas:57c3741f55d53490"
    
    tag "${meta.id}"

    publishDir "${params.outdir}/consensus_reports/${meta}"

    input:
    tuple val(meta), path(samplesheet), path(bed_files)
    val ucsc_session_url

    output:
    tuple val(meta), path("*"),           emit: results 

    script:
    """
    consensus.py \
        -s ${samplesheet} \
        -n ${meta.id} \
        -o consensus_results \
        -u "${ucsc_session_url}"
    """
    
    stub:
    """
    mkdir -p consensus_results
    touch consensus_results/${prefix}.txt
    """
}