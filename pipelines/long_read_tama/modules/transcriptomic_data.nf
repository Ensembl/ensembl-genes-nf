process TRANSCRIPTOMIC_DATA {
    tag "taxon ${taxon_id}:transcriptomic-data"
    label 'process_light'
    container { params.ensembl_genes_container ?: 'dockerhub.ebi.ac.uk/ensembl_genebuild/ensembl-genes-containers/ensembl-genes' }

    input:
    val taxon_id

    output:
    path 'transcriptomic_candidates.tsv', emit: candidates
    path 'versions.yml', emit: versions

    script:
    def tree = params.tree in [true, 'true'] ? '--tree' : ''
    """
    ${params.transcriptomic_data_command} -t ${taxon_id} -f transcriptomic_candidates.tsv -r long ${tree}
    printf '"%s":\\n    transcriptomic_data: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'stub\\tSRR000001\\t1\\tstub.fastq.gz\\t-1\\t1\\t0\\tENA\\tONT\\tstub\\thttps://example.org/SRR000001.fastq.gz\\td41d8cd98f00b204e9800998ecf8427e\\n' > transcriptomic_candidates.tsv
    printf '"%s":\\n    transcriptomic_data: stub\\n' '${task.process}' > versions.yml
    """
}
