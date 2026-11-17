process INDEX_AND_VALIDATE_BAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'
    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: indexed
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    genome=\$(find staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
    [[ -n "\$genome" && -s "\$genome" ]]
    samtools faidx "\$genome"
    for bam in staged/bam_*.bam; do
        [[ -s "\$bam" ]]
        samtools quickcheck -v "\$bam"
        samtools view -H "\$bam" | grep -q 'SO:coordinate' || { echo "BAM is not coordinate sorted: \$bam" >&2; exit 1; }
        bai="\${bam}.bai"
        if [ ! -s "\$bai" ]; then samtools index -@ ${task.cpus} "\$bam" "\$bai"; fi
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
"\${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    mkdir -p staged
    touch staged/bam_01.bam.bai staged/bam_02.bam.bai
    touch versions.yml
    """
}
