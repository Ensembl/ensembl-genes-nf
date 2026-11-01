process PARSE_GFFCOMPARE_TMAP {
    tag { "${meta.pair_id}_${direction}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/gffcompare" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), val(direction), path(tmap)

    output:
    tuple val(meta), val(direction), path("${direction}.classes.tsv"), emit: classes
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    python - '${tmap}' '${direction}.classes.tsv' '${meta.pair_id}' <<'PY'
    import csv, sys
    source, output, pair_id = sys.argv[1:]
    counts = {}
    total = 0
    with open(source, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\\t')
        for row in reader:
            code = row.get('class_code', '')
            counts[code] = counts.get(code, 0) + 1
            total += 1
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\\t', lineterminator='\\n')
        writer.writerow(['pair_id', 'direction', 'class_code', 'n_models', 'denominator', 'rate'])
        for code in sorted(counts):
            writer.writerow([pair_id, '${direction}', code, counts[code], total, f'{counts[code] / total:.6f}' if total else '0'])
    PY
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tdirection\\tclass_code\\tn_models\\tdenominator\\trate\\n' > '${direction}.classes.tsv'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
