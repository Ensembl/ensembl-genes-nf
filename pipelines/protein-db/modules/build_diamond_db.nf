process BUILD_DIAMOND_DB {
    label 'process_high'
    container params.use_containers ? params.tool_containers.diamond : null
    publishDir "${params.outdir}/databases", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${meta.id}/${filename}" }
    tag "${meta.id}:${meta.source}:${meta.scope}"

    input:
    tuple val(meta), path(protein_fasta), path(source_metadata)

    output:
    tuple val(meta), path('database.dmnd'), path('database_proteins.faa'), path('provenance.tsv'), emit: database
    path 'versions.yml', emit: versions

    script:
    """
    diamond makedb --in ${protein_fasta} --db database --threads ${task.cpus}
    cp ${protein_fasta} database_proteins.faa

    cp ${source_metadata} provenance.tsv
    printf 'diamond_version\\t%s\\n' "\$(diamond version | head -n1)" >> provenance.tsv
    printf 'diamond_database_sha256\\t%s\\n' "\$(sha256sum database.dmnd | awk '{print \\$1}')" >> provenance.tsv
    printf 'database_proteins_sha256\\t%s\\n' "\$(sha256sum database_proteins.faa | awk '{print \\$1}')" >> provenance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    """
    printf 'stub diamond database for ${meta.id}\\n' > database.dmnd
    cp ${protein_fasta} database_proteins.faa
    cp ${source_metadata} provenance.tsv
    printf 'diamond_version\\t2.1.24 (stub)\\n' >> provenance.tsv
    printf 'diamond_database_sha256\\tstub\\n' >> provenance.tsv
    printf 'database_proteins_sha256\\tstub\\n' >> provenance.tsv
    printf '"%s":\\n    diamond: 2.1.24\\n' '${task.process}' > versions.yml
    """
}
