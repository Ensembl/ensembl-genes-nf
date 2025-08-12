
process FINAL_VERIFICATION_DATACHECKS {

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/final_verification",
        pattern: "*.{txt,tap,log}"

    input:
    tuple val(meta), val(core_name)
    
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
        -output_file "\$TAPFILE" >> "\$LOGFILE" 2>&1

    if [ -s "\$TAPFILE" ]; then
        # Extract overall PASS/FAIL summary from TAP
        grep -E "ok |not ok" "\$TAPFILE" > "\$SUMMARY" || echo "No results found" > "\$SUMMARY"
    else
        echo "No TAP file produced" > "\$SUMMARY"
    fi
    """
}
