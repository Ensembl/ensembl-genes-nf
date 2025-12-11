// TODO: Rename this process to match your tool/module name (e.g., MY_TOOL)
process REPORT_CONSENSUS {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pip_duckdb_pandas:57c3741f55d53490"
    
    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/consensus_reports/${meta}"

    input:
    tuple val(meta), path(samplesheet), path(bed_files)
    val ucsc_session_url

    output:
    tuple val(meta), path("consensus_results/"),             emit: results 


    script:
    """
    consensus.py \
        -s ${samplesheet} \
        -n ${meta} \
        -o consensus_results \
        -u ${ucsc_session_url}
    """
    
    // TODO: Update stub section to create mock output files matching your tool's actual outputs
    stub:
    """
    mkdir -p consensus_results
    touch consensus_results/${prefix}.txt
    """
}