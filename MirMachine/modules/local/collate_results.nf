

process COLLATE_RESULTS {
    tag "Preprocessing miRNA data files"
    
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    path(input_files, stageAs: 'heatmap_*/*.heatmap.csv')

    output:
    path "${date}_miRNA_heatmap_for_R.csv", emit: heatmap
    path "${date}_miRNA_metadata.csv", emit: metadata
    
    script:
    date = new java.util.Date().format('yyyyMMdd')
    """
    # Collate all input files and extract required columns
    find heatmap_* -name '*.heatmap.csv' | xargs cat | cut -f1,2,7,8 | grep -v species | awk 'BEGIN{print "species,mode,family,node,tgff,filtered"} {sub(".PRE",""); print}' > ${date}_miRNA_heatmap.csv
    
    # Create R input file by removing comment lines
    sed '/^#/d' ${date}_miRNA_heatmap.csv > ${date}_miRNA_heatmap_for_R.csv
    
    # Create metadata file with species and family count info
    grep -E '# Total families searched:|# Species:' ${date}_miRNA_heatmap.csv | awk 'NR%2{printf "%s,",\$0;next;}1' | tr ' ' ',' > ${date}_miRNA_metadata.csv
    """
}