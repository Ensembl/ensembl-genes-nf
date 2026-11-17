process VALIDATE_GENOME {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_light'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: validated

    script:
    """
    genome=\$(find staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
    [[ -n "\$genome" && -s "\$genome" ]]
    if [[ "\$genome" == *.gz ]]; then
        gzip -cd "\$genome"
    else
        cat "\$genome"
    fi | awk '!/^>/ && /[a-z]/ { found=1; exit } END { if (!found) exit 1 }' || {
        echo "Genome appears not to be soft-masked: no lowercase bases found" >&2
        exit 1
    }
    """

    stub:
    """
    mkdir -p staged
    """
}
