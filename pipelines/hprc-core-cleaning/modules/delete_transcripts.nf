
/*
 * Problem: Specific transcripts identified as broken (missing translations, invalid 
 * coordinates, or other structural issues) need to be removed from the database.
 *
 * Solution: Use Ensembl's delete_transcripts.pl script to safely remove problematic 
 * transcripts while maintaining referential integrity and cascading deletions.
 */

process DELETE_TRANSCRIPTS {

    tag "${meta.id}"

    // errorStrategy 'ignore'

    publishDir "${params.outdir}/delete_transcripts",
        pattern: "*.{txt,log}"

    input:
    tuple val(meta), val(core_name), path(file)
    val(db_pass)

    output:
    tuple val(meta), path("${core_name}.txt"), emit: results
    tuple val(meta), path("${core_name}.log"), emit: logs

    script:
    """
    LOGFILE="${core_name}.log"
    SUMMARY="${core_name}.txt"

    echo "=== Deleting transcripts for ${meta.id} (${core_name}) ===" > "\$LOGFILE"
    echo "Core DB: ${core_name}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    echo "Transcript list: \$(wc -l < "$file") lines" >> "\$LOGFILE"

    if [ "\$(wc -l < "$file")" -eq 0 ]; then
        echo "No transcripts to delete." >> "\$LOGFILE"
        echo "No deletions" > "\$SUMMARY"
        exit 0
    elif [ "\$(wc -l < "$file")" -gt 10 ]; then
        echo "Too many transcripts to delete. Suspect bigger issue." >> "\$LOGFILE"
        echo "Too many transcripts" > "\$SUMMARY"
        exit 1
    fi

    perl /hps/software/users/ensembl/genebuild/jackt/modenv/hprc2-prod/ensembl-analysis/scripts/genebuild/delete_transcripts.pl \
        --dbhost mysql-ens-genebuild-prod-1 \
        --dbport 4527 --dbuser ensadmin --dbpass "$db_pass" \
        --dbname "$core_name" \
        --stable_id "$file" >> "\$LOGFILE" 2>&1

    echo "Deletions complete" > "\$SUMMARY"
    """
}
