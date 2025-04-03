

process GENERATE_SCORES {
    tag "Calculating microRNA scores"

    publishDir "${params.outdir}/mirmachine/", mode: 'copy'
    // errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    input:
    path heatmap_file
    path metadata_file
    
    output:
    path "${date}_microRNA_score_output.tsv", emit: scores
    
    script:
    date = new java.util.Date().format('yyyyMMdd')
    """
    microRNA_scoring.py --input ${heatmap_file} --metadata ${metadata_file} --output "${date}_microRNA_score_output.tsv"
    """
}