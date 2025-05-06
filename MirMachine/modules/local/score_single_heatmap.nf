


process GENERATE_SINGLE_SCORES {

    publishDir "${params.outdir}/mirmachine/", mode: 'copy'
    // errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    input:
    path heatmap_file

    
    output:
    path "${heatmap_file.basename}_scores.tsv ", emit: scores
    
    script:
    """
    python single_heatmap_miRNA_scoring.py -i $heatmap_file -o ${heatmap_file.basename}_scores.tsv  

    """
}