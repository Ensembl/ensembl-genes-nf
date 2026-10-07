process SAMTOOLS_QUICKCHECK {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes)

    output:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes), emit: checked
    path "${meta.safe_id}.bam_qc.versions.yml", emit: versions

    script:
    """
    echo 'Checking BAM readability and coordinate sort for ${meta.gca_accession}:${meta.id}' >&2
    samtools quickcheck -v ${bam}
    samtools view -H ${bam} | grep -q 'SO:coordinate' || { echo 'BAM is not coordinate sorted' >&2; exit 1; }
    test -s ${bai} || { echo 'BAM index is empty' >&2; exit 1; }
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > ${meta.safe_id}.bam_qc.versions.yml
    """

    stub:
    """
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > ${meta.safe_id}.bam_qc.versions.yml
    """
}
