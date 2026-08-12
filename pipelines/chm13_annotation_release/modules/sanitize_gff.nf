process SANITIZE_GFF {
    tag "${meta.id}"
    label 'process_low'

    conda 'conda-forge::python=3.12'

    input:
    tuple val(meta), path(projected_gff), path(manual_gff), path(decision_tsv), path(assembly_fasta), path(assembly_fai)
    path sanitize_script

    output:
    tuple val(meta), path("${meta.id}.projected.sanitized.gff3"), path("${meta.id}.manual.sanitized.gff3"), path(decision_tsv), path(assembly_fasta), path(assembly_fai), emit: sources
    path "${meta.id}.gff_sanitization.*.json", emit: reports

    script:
    """
    set -euo pipefail
    python3 ${sanitize_script} --input ${projected_gff} --output ${meta.id}.projected.sanitized.gff3 --kind projected --report ${meta.id}.gff_sanitization.projected.json
    python3 ${sanitize_script} --input ${manual_gff} --output ${meta.id}.manual.sanitized.gff3 --kind manual --report ${meta.id}.gff_sanitization.manual.json
    """

    stub:
    """
    touch ${meta.id}.projected.sanitized.gff3
    touch ${meta.id}.manual.sanitized.gff3
    touch ${meta.id}.gff_sanitization.projected.json
    touch ${meta.id}.gff_sanitization.manual.json
    """
}
