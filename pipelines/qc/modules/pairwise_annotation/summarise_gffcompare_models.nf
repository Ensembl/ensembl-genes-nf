process SUMMARISE_GFFCOMPARE_MODELS {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/gffcompare" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(classes)

    output:
    tuple val(meta), path('gffcompare_summary.tsv'), emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    python - ${classes} '${meta.pair_id}' > gffcompare_summary.tsv <<'PY'
    import csv, sys
    writer = csv.writer(sys.stdout, delimiter='\\t', lineterminator='\\n')
    writer.writerow(['pair_id', 'direction', 'n_models', 'n_class_codes'])
    for path in sys.argv[1:-1]:
        with open(path, newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\\t'))
        direction = rows[0]['direction'] if rows else path
        writer.writerow([sys.argv[3], direction, sum(int(row.get('n_models', 0)) for row in rows), len(rows)])
    PY
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tdirection\\tn_models\\tn_class_codes\\n' > gffcompare_summary.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
