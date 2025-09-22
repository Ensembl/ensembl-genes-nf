
/*
 * Problem: DNA sequences, attributes, and synonyms reference deleted seq_regions, 
 * creating orphaned records that violate referential integrity.
 *
 * Solution: Delete all orphaned DNA, seq_region_attrib, and seq_region_synonym 
 * records that point to non-existent seq_region entries.
 */

process TIDY_DUPLICATE_COORD_SYSTEMS {

    tag "${meta.id}"

    // errorStrategy 'ignore'

    publishDir "${params.outdir}/tidy_duplicate_coord_systems",
        pattern: "*.{txt,log}"

    input:
    tuple val(meta), val(core_db)

    output:
    tuple val(meta), path("${core_db}.metadata.log"), emit: logs

    script:
    """
    LOGFILE="${core_db}.metadata.log"
    SUMMARY="${core_db}.metadata.txt"

    echo "=== Updating metadata for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Server: ${meta.server}" >> "\$LOGFILE"
    echo "Core DB: ${core_db}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    # Run queries individually, log result
    ${meta.server} ${core_db} -e "DELETE FROM dna where seq_region_id not in (SELECT seq_region_id FROM seq_region);"
    ${meta.server} ${core_db} -e "DELETE FROM seq_region_attrib where seq_region_id not in (SELECT seq_region_id FROM seq_region);"
    ${meta.server} ${core_db} -e "DELETE FROM seq_region_synonym where seq_region_id not in (SELECT seq_region_id FROM seq_region);"

    # Short summary: just count successful updates
    grep "Query OK" "\$LOGFILE" > "\$SUMMARY" || echo "No successful updates" > "\$SUMMARY"
    """
}
