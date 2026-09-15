process DISCOVER_LONG_READ_DATA {
    tag "taxon ${taxon_id}"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11--he2b4eab_0'

    input:
    val taxon_id
    path discoverer
    val cache_dir

    output:
    path 'discovery', emit: results

    script:
    def tree = params.tree ? '--tree' : ''
    def probe = params.discovery_probe ? '' : '--no-probe'
    """
    mkdir -p discovery
    python3 ${discoverer} ${taxon_id} discovery --cache-dir ${cache_dir} --target ${params.target_sample_count} --soft-download-budget ${params.soft_download_budget} --soft-raw-subread-budget ${params.soft_raw_subread_budget} ${tree} ${probe}
    """

    stub:
    """
    mkdir -p discovery
    printf 'discovered_runs\\t1\\nbiological_samples\\t1\\nselected_runs\\t1\\n' > discovery/selection_summary.json
    printf 'run_accession\\tselection_status\\nSRR000001\\tSELECTED\\n' > discovery/long_read_inventory.tsv
    printf 'run_accession\\tselection_status\\nSRR000001\\tSELECTED\\n' > discovery/proposed_long_read_manifest.tsv
    """
}
