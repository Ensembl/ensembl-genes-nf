process VALIDATE_SPLICE_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:482--hdc0a859_0'

    input:
    tuple val(meta), path(bigbed), path(bam)

    output:
    tuple val(meta), path(bigbed), path(bam), emit: validated
    path "${meta.safe_id}.splice_bigbed_validation.tsv", emit: validation

    script:
    """
    test -s ${bigbed} || { echo 'Splice-junction BigBed validation failed' >&2; exit 1; }
    bigBedInfo -as ${bigbed} > ${meta.safe_id}.embedded_splice_schema.as
    for field in chrom chromStart chromEnd name score strand unique_reads multi_reads splice_motif; do
        grep -q "[[:space:]]${field};" ${meta.safe_id}.embedded_splice_schema.as || {
            echo "Splice BigBed is missing embedded AutoSql field: ${field}" >&2
            exit 1
        }
    done
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tsplice_junction\\tcomplete\\n' > ${meta.safe_id}.splice_bigbed_validation.tsv
    """

    stub:
    """
    cat > ${meta.safe_id}.embedded_splice_schema.as <<'AUTOSQL'
    table spliceJunction
    "Splice-junction intervals and read support. Coordinates use 0-based, half-open BED intervals."
    (
        string chrom; "Reference sequence name"
        uint chromStart; "Junction start, 0-based"
        uint chromEnd; "Junction end, exclusive"
        string name; "Junction identifier"
        uint score; "Display score, 0-1000; capped from unique_reads"
        char[1] strand; "+, -, or ."
        uint unique_reads; "Number of uniquely mapping reads spanning the junction"
        uint multi_reads; "Number of multimapping reads spanning the junction"
        string splice_motif; "Splice-site dinucleotide pair; . means unavailable or non-canonical in STAR"
    )
    AUTOSQL
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tsplice_junction\\tcomplete\\n' > ${meta.safe_id}.splice_bigbed_validation.tsv
    """
}
