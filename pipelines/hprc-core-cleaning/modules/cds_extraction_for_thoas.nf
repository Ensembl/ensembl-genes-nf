process EXTRACT_CDS {

    tag "${core_name}"

    publishDir "${params.outdir}/extract_cds", pattern: "*"

    input:
    tuple val(meta), val(core_name)

    output:
    tuple val(meta), path("${core_name}_cds.fa"), emit: cds
    tuple val(meta), path("${core_name}.log"), emit: logs

    script:
    """
    LOGFILE="${core_name}.log"
    OUTFILE="${core_name}_cds.fa"

    echo "=== Extracting CDS for ${meta.id} (${core_name}) ===" > "\$LOGFILE"
    echo "Server: ${meta.server}" >> "\$LOGFILE"
    echo "Core DB: ${core_name}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    perl /hps/software/users/ensembl/applications/bilal/thoas/ensembl-core-mongodb-loading/src/ensembl/extract_cds_from_ens.pl \
        --host=mysql-ens-genebuild-prod-1.ebi.ac.uk \
        --user=ensro \
        --port=4527 \
        --species=homo_sapiens \
        --assembly='${core_name}' \
        --database="${core_name}" \
        > "\$OUTFILE" 2>> "\$LOGFILE"

    COUNT=\$(grep -c '^>' "\$OUTFILE" || true)
    echo "Total CDS sequences extracted: \$COUNT" >> "\$LOGFILE"
    """
}
