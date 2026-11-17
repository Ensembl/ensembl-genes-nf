process MERGE_BAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: merged
    path 'versions.yml', emit: versions

    script:
    """
    bam_inputs=(staged/bam_*.bam)
    [[ -e "\${bam_inputs[0]}" ]]
    samtools merge -@ ${task.cpus} -o staged/all_runs.merged.bam "\${bam_inputs[@]}"
    samtools sort -@ ${task.cpus} -o staged/all_runs.merged.sorted.bam staged/all_runs.merged.bam
    samtools index -@ ${task.cpus} staged/all_runs.merged.sorted.bam
    cat > versions.yml <<VERSIONS
"\${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    touch staged/all_runs.merged.bam staged/all_runs.merged.sorted.bam staged/all_runs.merged.sorted.bam.bai
    touch versions.yml
    """
}
