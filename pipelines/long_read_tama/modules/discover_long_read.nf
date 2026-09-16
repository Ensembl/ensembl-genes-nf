process DISCOVER_LONG_READ_DATA {
    tag "taxon ${taxon_id}"
    label 'process_light'
    container { params.ensembl_genes_container ?: 'dockerhub.ebi.ac.uk/ensembl_genebuild/ensembl-genes-containers/ensembl-genes' }

    input:
    val taxon_id
    path discoverer
    path classifier
    path resolver
    val cache_dir

    output:
    path 'discovery', emit: results
    path 'versions.yml', emit: versions

    script:
    def tree = params.tree ? '--tree' : ''
    def probe = params.discovery_probe ? '' : '--no-probe'
    """
    mkdir -p discovery
    ${params.transcriptomic_data_command} -t ${taxon_id} -f discovery/transcriptomic_candidates.tsv -r long ${tree}
    python3 ${discoverer} ${taxon_id} discovery --cache-dir ${cache_dir} --candidate-file discovery/transcriptomic_candidates.tsv --target ${params.target_sample_count} --soft-download-budget ${params.soft_download_budget} --soft-raw-subread-budget ${params.soft_raw_subread_budget} ${tree} ${probe}
    printf '"%s":\\n    python: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p discovery
    printf 'stub\tSRR000001\t1\tstub.fastq.gz\t-1\t1\t0\tENA\tONT\tstub\thttps://example.org/SRR000001.fastq.gz\td41d8cd98f00b204e9800998ecf8427e\n' > discovery/transcriptomic_candidates.tsv
    printf 'discovered_runs\\t1\\nbiological_samples\\t1\\nselected_runs\\t1\\n' > discovery/selection_summary.json
    printf 'run_accession\\tselection_status\\nSRR000001\\tSELECTED\\n' > discovery/long_read_inventory.tsv
    printf 'run_accession\\ttissue\\tfilename\\turl\\tmd5\\tsource\\tplatform\\tdescription\\tselection_status\\tselection_reason\\nSRR000001\\tstub\\tSRR000001.fastq.gz\\thttps://example.org/SRR000001.fastq.gz\\td41d8cd98f00b204e9800998ecf8427e\\tENA\\tOXFORD_NANOPORE\\tstub\\tSELECTED\\tstub\\n' > discovery/proposed_long_read_manifest.tsv
    printf '"%s":\\n    discovery: stub\\n' '${task.process}' > versions.yml
    """
}
