def shell_quote(value) {
    "'${value.toString().replace("'", "'\\\"'\\\"'")}'"
}

process FETCH_UNIPROT_FASTA {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${meta.source}:${meta.scope}"

    input:
    tuple val(meta), val(build)
    val uniprot_api
    val page_size
    path fetch_script

    output:
    tuple val(meta), path('proteins.faa'), path('source_metadata.tsv'), emit: fasta
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${fetch_script} \\
        --api ${shell_quote(uniprot_api)} \\
        --build-id ${shell_quote(build.build_id)} \\
        --source ${shell_quote(build.source)} \\
        --scope ${shell_quote(build.scope)} \\
        --taxon-id ${shell_quote(build.taxon_id)} \\
        --taxon-name ${shell_quote(build.taxon_name)} \\
        --query-extra ${shell_quote(build.query_extra ?: '')} \\
        --include-isoforms ${build.include_isoforms ? 'true' : 'false'} \\
        --page-size ${page_size} \\
        --fasta proteins.faa \\
        --metadata source_metadata.tsv

    printf '"%s":\\n    python: %s\\n' '${task.process}' "\$(python3 --version 2>&1 | sed 's/^Python //')" > versions.yml
    """

    stub:
    """
    printf '>stub|${meta.id}|STUB_PROTEIN_${meta.source}\\nMPEPTIDE\\n' > proteins.faa
    printf 'build_id\\t%s\\nsource\\t%s\\nscope\\t%s\\ntaxon_id\\t%s\\ntaxon_name\\t%s\\nquery\\tstub\\nsequence_count\\t1\\nfasta_sha256\\tstub\\n' \\
        '${build.build_id}' '${build.source}' '${build.scope}' '${build.taxon_id}' '${build.taxon_name}' > source_metadata.tsv
    printf '"%s":\\n    python: 3.11\\n' '${task.process}' > versions.yml
    """
}
