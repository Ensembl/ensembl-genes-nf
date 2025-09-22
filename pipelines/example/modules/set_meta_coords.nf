
/*
 * Problem: The meta_coord table contains incorrect max_length values for genomic 
 * features, likely from template defaults that don't match the actual data.
 *
 * Solution: Calculate actual maximum feature lengths from transcript, gene, exon, 
 * and repeat_feature tables and update meta_coord entries accordingly.
 */

process SET_META_CORDS {

    tag "${meta.id}"

    // errorStrategy 'ignore'

    publishDir "${params.outdir}/metadcoords_updates"

    input:
    tuple val(meta), val(core_db)

    output:
    tuple val(meta), path("${core_db}.metacoord.txt"), emit: results
    tuple val(meta), path("${core_db}.metacoord.log"), emit: logs

    script:
    """
    LOGFILE="${core_db}.metacoord.log"
    SUMMARY="${core_db}.metacoord.txt"

    echo "=== Updating meta coords for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Server: ${meta.server}" >> "\$LOGFILE"
    echo "Core DB: ${core_db}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    # First, check current state
    echo "Current meta_coord state:" >> "\$LOGFILE"
    ${meta.server} ${core_db} -e "SELECT table_name, max_length FROM meta_coord WHERE table_name = 'transcript';" >> "\$LOGFILE" 2>&1
    
    # Calculate the max transcript length
    echo "Calculating max transcript length..." >> "\$LOGFILE"
    tx_MAX_LENGTH=\$(${meta.server} ${core_db} -e "SELECT MAX(seq_region_end - seq_region_start + 1) FROM transcript;" 2>> "\$LOGFILE" | tail -n 1)
    gene_MAX_LENGTH=\$(${meta.server} ${core_db} -e "SELECT MAX(seq_region_end - seq_region_start + 1) FROM gene;" 2>> "\$LOGFILE" | tail -n 1)
    exon_MAX_LENGTH=\$(${meta.server} ${core_db} -e "SELECT MAX(seq_region_end - seq_region_start + 1) FROM exon;" 2>> "\$LOGFILE" | tail -n 1)
    repeat_MAX_LENGTH=\$(${meta.server} ${core_db} -e 'SELECT MAX(seq_region_end - seq_region_start + 1) FROM `repeat_feature`;' 2>> "\$LOGFILE" | tail -n 1)
    
    if [[ "\$tx_MAX_LENGTH" =~ ^[0-9]+\$ ]]; then
        echo "Max transcript length: \$tx_MAX_LENGTH" >> "\$LOGFILE"
        echo "Max gene length: \$gene_MAX_LENGTH" >> "\$LOGFILE"
        echo "Max exon length: \$exon_MAX_LENGTH" >> "\$LOGFILE"
        echo "Max repeat length: \$repeat_MAX_LENGTH" >> "\$LOGFILE"
        
        # Update meta_coord table
        ${meta.server} ${core_db} -e "UPDATE meta_coord SET max_length = \$tx_MAX_LENGTH WHERE table_name = 'transcript';" >> "\$LOGFILE" 2>&1
        ${meta.server} ${core_db} -e "UPDATE meta_coord SET max_length = \$gene_MAX_LENGTH WHERE table_name = 'gene';" >> "\$LOGFILE" 2>&1
        ${meta.server} ${core_db} -e "UPDATE meta_coord SET max_length = \$exon_MAX_LENGTH WHERE table_name = 'exon';" >> "\$LOGFILE" 2>&1
        ${meta.server} ${core_db} -e "UPDATE meta_coord SET max_length = \$repeat_MAX_LENGTH WHERE table_name = 'repeat_feature';" >> "\$LOGFILE" 2>&1
        # Remove any existing entries for dna_align_feature
        ${meta.server} ${core_db} -e "DELETE FROM meta_coord WHERE table_name = 'dna_align_feature';" >> "\$LOGFILE" 2>&1

        # Verify the update
        echo "Updated meta_coord state:" >> "\$LOGFILE"
        ${meta.server} ${core_db} -e "SELECT table_name, max_length FROM meta_coord WHERE table_name = 'transcript';" >> "\$LOGFILE" 2>&1
        
        echo "Successfully updated transcript max_length to \$tx_MAX_LENGTH" > "\$SUMMARY"
    else
        echo "ERROR: Could not calculate max transcript length" >> "\$LOGFILE"
        echo "ERROR: Failed to calculate max transcript length" > "\$SUMMARY"
    fi
    """
}

