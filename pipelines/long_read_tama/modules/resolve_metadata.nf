process RESOLVE_LONG_READ_METADATA {
    tag 'ENA metadata inventory'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    path manifest
    path resolver
    val cache_dir

    output:
    path 'metadata.json', emit: metadata
    path 'metadata_cache', emit: cache

    script:
    """
    mkdir -p ${cache_dir}
    tail -n +2 ${manifest} | cut -f1 > accessions.txt
    python3 ${resolver} accessions.txt ${cache_dir} metadata.json
    """

    stub:
    """
    mkdir -p metadata_cache
    printf '{"SRR000001": {"instrument_platform": "PACBIO_SMRT", "instrument_model": "stub", "sra_platform": "unavailable", "sra_spot_group": "unavailable", "submitted_ftp": "", "submitted_md5": "", "submitted_format": "", "fastq_ftp": "https://example.org/SRR000001.fastq.gz", "fastq_md5": "d41d8cd98f00b204e9800998ecf8427e"}}\\n' > metadata.json
    """
}
