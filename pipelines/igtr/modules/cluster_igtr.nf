// CLUSTER_IGTR
// Merge GFF3 files from all GenBlast batches, cluster overlapping gene models
// on the same chromosome+strand, and retain the best model per locus
// (highest combined PID + Coverage score). Emulates HiveCollapseIGTR.

process CLUSTER_IGTR {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/igtr", mode: 'copy', pattern: "*.clustered.gff3"

    input:
    tuple val(meta), path(gff3_files)   // all filtered GFF3s collected

    output:
    tuple val(meta), path("*.clustered.gff3"), emit: gff3
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    # Merge all input GFF3 files (skip per-file headers)
    echo '##gff-version 3' > merged.gff3
    cat ${gff3_files} | grep -v '^#' >> merged.gff3

    cluster_igtr.py \\
        --gff3      merged.gff3 \\
        --out       ${prefix}.clustered.gff3 \\
        --sample-id ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    cat > ${prefix}.clustered.gff3 << 'STUB_GFF3'
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t85.5\t+\t.\tID=${prefix}_igtr_gene_00000001;Name=STUB;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=${prefix}_igtr_transcript_00000001;Parent=${prefix}_igtr_gene_00000001;Name=STUB;PID=85.50;Coverage=98.30;biotype=ig_gene
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=${prefix}_igtr_exon_00000001;Parent=${prefix}_igtr_transcript_00000001
    STUB_GFF3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
