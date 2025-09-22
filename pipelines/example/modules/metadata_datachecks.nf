
/*
 * Problem: After database cleanup and metadata updates, need to verify that all 
 * critical integrity checks pass (foreign keys, display xrefs, meta coordinates).
 *
 * Solution: Run Ensembl datachecks on key validation tests to confirm the database 
 * is properly structured and ready for production use.
 */

process FINAL_VERIFICATION_DATACHECKS {

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/final_verification",
        pattern: "*.{txt,tap,log}"

    input:
    tuple val(meta), val(core_name)
    tuple val(meta), path(metadata_results)
    tuple val(meta), path(metacoords_results)
    
    output:
    tuple val(meta), path("${core_name}.txt"), emit: results
    tuple val(meta), path("${core_name}.tap"), emit: tap
    tuple val(meta), path("${core_name}.log"), emit: logs

    script:
    """
    LOGFILE="${core_name}.log"
    SUMMARY="${core_name}.txt"
    TAPFILE="${core_name}.tap"

    echo "=== Running final verification datachecks for ${meta.id} (${core_name}) ===" > "\$LOGFILE"
    echo "Core DB: ${core_name}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    run_datachecks.pl \\
        -H mysql-ens-genebuild-prod-1 \\
        -P 4527 \\
        -u ensro \\
        -dbname "${core_name}" \\
        -names DisplayXrefExists,ForeignKeys,MetaCoord,MetaKeyCardinality,MetaKeyConditional,MetaKeyFormat \\
        -output_file "\$TAPFILE" >> "\$LOGFILE" 

    if [ -s "\$TAPFILE" ]; then
        # Extract overall PASS/FAIL summary from TAP
        grep -E "ok |not ok" "\$TAPFILE" > "\$SUMMARY" || echo "No results found" > "\$SUMMARY"
    else
        echo "No TAP file produced" > "\$SUMMARY"
    fi
    """
}
