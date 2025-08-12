

process FIND_MISSING_TRANSLATIONS {

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/missing_translations"

    input:
    tuple val(meta), val(core_db)

    output:
    tuple val(meta), val(core_db), path("${core_db}.txt"), emit: results
    tuple val(meta), path("${core_db}.log"), emit: logs

    script:
    """
    LOGFILE="${core_db}.log"
    
    echo "=== Checking protein coding transcripts for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Server: ${meta.server}" >> "\$LOGFILE"
    echo "Core DB: ${core_db}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    # Run queries individually, log result
    ${meta.server} ${core_db} -NB -e "SELECT stable_id FROM transcript WHERE biotype = 'protein_coding' AND transcript_id NOT IN (SELECT transcript_id FROM translation);"  > ${core_db}.txt 
    if [ -s ${core_db}.txt ]; then
        echo "Found missing translations for the following transcripts:" >> "\$LOGFILE"
        cat ${core_db}.txt >> "\$LOGFILE"
    else
        echo "No missing translations found." >> "\$LOGFILE"
    fi
    """
}
