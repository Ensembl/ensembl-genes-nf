

/*
 * Problem: Translation records have seq_start/seq_end coordinates that exceed 
 * the length of their corresponding start/end exons, creating invalid CDS boundaries.
 * These fixes are required to esure thoas loading.
 *
 * Solution: Identify translations where coordinates are out-of-bounds relative to 
 * exon lengths and report transcript stable_ids for coordinate correction.
 */


process FIND_BROKEN_TRANSLATIONS {

    tag "${meta.id}"

    // errorStrategy 'ignore'

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
    ${meta.server} ${core_db} -NB -e "SELECT transcript.stable_id FROM translation t JOIN transcript ON t.transcript_id = transcript.transcript_id JOIN exon ex_start ON t.start_exon_id = ex_start.exon_id JOIN exon ex_end ON t.end_exon_id = ex_end.exon_id WHERE t.seq_start > (ex_start.seq_region_end - ex_start.seq_region_start + 1) OR t.seq_end > (ex_end.seq_region_end - ex_end.seq_region_start + 1);"  > ${core_db}.txt 
    if [ -s ${core_db}.txt ]; then
        echo "Found missing translations for the following transcripts:" >> "\$LOGFILE"
        cat ${core_db}.txt >> "\$LOGFILE"
    else
        echo "No missing translations found." >> "\$LOGFILE"
    fi
    """
}
