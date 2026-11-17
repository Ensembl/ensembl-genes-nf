process VALIDATE_BAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: validated

    script:
    """
    for bam in staged/bam_*.bam; do
        [[ -s "\$bam" ]]
        samtools quickcheck -v "\$bam"
        samtools view -H "\$bam" | grep -q 'SO:coordinate' || {
            echo "BAM is not coordinate sorted: \$bam" >&2
            exit 1
        }
    done
    """

    stub:
    """
    mkdir -p staged
    """
}
