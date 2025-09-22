
/*
 * Problem: Database metadata contains placeholder values (YYYY-MM, '0') and 
 * incorrect/outdated entries from genome projection templates.
 *
 * Solution: Update metadata fields with correct dates, build IDs, provider 
 * information, and remove obsolete strain/repeat analysis entries.
 */

process SET_MISSING_METADATA {

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/metadata_updates"

    input:
    tuple val(meta), val(core_db)

    output:
    tuple val(meta), path("${core_db}.metadata.txt"), emit: results
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
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = '2025-08' WHERE meta_key = 'genebuild.initial_release_date' AND meta_value = 'YYYY-MM';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = '2025-08' WHERE meta_key = 'genebuild.last_geneset_update' AND meta_value = 'YYYY-MM';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = '61' WHERE meta_key = 'genebuild.id' AND meta_value = '0';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = 'HPRC' WHERE meta_key = 'assembly.provider_name';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = 'https://humanpangenome.org' WHERE meta_key = 'assembly.provider_url';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = 'Ensembl' WHERE meta_key = 'genebuild.provider_name';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = 'https://beta.ensembl.org/help/articles/human-genome-automated-annotation' WHERE meta_key = 'genebuild.provider_url';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = '' WHERE meta_key = 'species.strain' AND meta_value = 'reference';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE meta SET meta_value = '' WHERE meta_key = 'strain.type' AND meta_value = 'strain';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "DELETE FROM meta WHERE meta_key = 'strain.type';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "DELETE FROM meta WHERE meta_key = 'species.strain';"  >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "DELETE FROM meta WHERE meta_key = 'repeat.analysis' AND meta_value = 'repeatdetector';"  >> "\$LOGFILE" 2>&1

    # Short summary: just count successful updates
    grep "Query OK" "\$LOGFILE" > "\$SUMMARY" || echo "No successful updates" > "\$SUMMARY"
    """
}
