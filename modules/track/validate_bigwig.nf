process VALIDATE_BIGWIG {
    tag { "${meta.gca_accession}:${meta.id}:${meta.coverage_strand}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bigwiginfo:482--h0b57e2e_0'

    input:
    tuple val(meta), path(bigwig)

    output:
    tuple val(meta), path(bigwig), emit: validated
    path "${meta.safe_id}.${meta.coverage_strand}.bigwig_validation.tsv", emit: validation
    path "${meta.safe_id}.${meta.coverage_strand}.bigwig_info.txt", emit: info
    path "${meta.safe_id}.${meta.coverage_strand}.bigwig_validation.versions.yml", emit: versions

    script:
    """
    test -s ${bigwig} || { echo 'BigWig validation failed: empty file' >&2; exit 1; }
    bigWigInfo ${bigwig} > ${meta.safe_id}.${meta.coverage_strand}.bigwig_info.txt
    grep -q '^version' ${meta.safe_id}.${meta.coverage_strand}.bigwig_info.txt || {
        echo 'BigWig validation failed: bigWigInfo did not recognize the file' >&2
        exit 1
    }
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tcoverage_${meta.coverage_strand}\\tcomplete\\n' > ${meta.safe_id}.${meta.coverage_strand}.bigwig_validation.tsv
    printf '"%s":\\n    bigWig_validation: repository-shell\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.bigwig_validation.versions.yml
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tcoverage_${meta.coverage_strand}\\tcomplete\\n' > ${meta.safe_id}.${meta.coverage_strand}.bigwig_validation.tsv
    printf 'stub BigWig information\\nversion 482\\n' > ${meta.safe_id}.${meta.coverage_strand}.bigwig_info.txt
    printf '"%s":\\n    bigWig_validation: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.bigwig_validation.versions.yml
    """
}
