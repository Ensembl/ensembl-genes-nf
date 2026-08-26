process PREPARE_RPBP_GENOME {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0'
    input:
    tuple val(meta), path(fastq)
    path gtf
    path fasta
    path ribosomal_fasta
    path adapter_fasta
    output:
    tuple val(meta), path('raw/rpbp.yaml'), emit: config
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    mkdir -p raw
    cat > raw/rpbp.yaml <<-END_CONFIG
    gtf: ${gtf}
    fasta: ${fasta}
    ribosomal_fasta: ${ribosomal_fasta}
    genome_name: ${meta.id}_rpbp
    genome_base_path: raw/genome
    ribosomal_index: raw/ribosomal_index
    star_index: raw/star_index
    riboseq_samples:
      ${meta.id}: ${fastq}
    adapter_file: ${adapter_fasta}
    riboseq_data: raw/data
    END_CONFIG
    prepare-rpbp-genome raw/rpbp.yaml --star-options '--genomeSAindexNbases 8' --num-cpus ${task.cpus ?: 1} --logging-level INFO
    printf '"%s":\n    Rp-Bp: 3.0.1\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf 'stub\n' > raw/rpbp.yaml
    printf '"stub":\n    Rp-Bp: stub\n' > versions.yml
    """
}
