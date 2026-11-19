process INDEX_AND_VALIDATE_BAMS {
    tag "${meta.assembly_accession}:${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: indexed
    path 'versions.yml', emit: versions

    script:
    """
    genome=\$(find -L staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
    [[ -n "\$genome" && -s "\$genome" ]]
    samtools faidx "\$genome"
    for cram in staged/cram_*.cram; do
        [[ -s "\$cram" ]]
        samtools quickcheck -v "\$cram"
    done
    if [[ "\$genome" == *.gz ]]; then
        gzip -cd "\$genome"
    else
        cat "\$genome"
    fi | awk '!/^>/ && /[a-z]/ { found=1; exit } END { if (!found) exit 1 }' || {
        echo "Genome appears not to be soft-masked: no lowercase bases found" >&2
        exit 1
    }
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    mkdir -p staged
    touch staged/cram_01.cram staged/cram_02.cram
    touch versions.yml
    """
}
