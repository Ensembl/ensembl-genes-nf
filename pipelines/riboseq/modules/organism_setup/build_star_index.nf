/*
 * BUILD_STAR_INDEX
 * Build STAR genome index for alignment
 */

process BUILD_STAR_INDEX {
    tag "${organism}_${version}"
    label 'process_high'

    container "quay.io/biocontainers/star:2.7.6a--0"

    publishDir "${params.outdir}/organism_setup", mode: 'copy'

    input:
    path(genome_fasta)
    path(gtf)
    val(organism)
    val(version)

    output:
    path "star_index",     emit: index
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Only set memory limit if explicitly provided and reasonable (>= 8GB)
    def memory = (task.memory && task.memory.toGiga() >= 8) ?
        "--limitGenomeGenerateRAM ${task.memory.toBytes()}" : ''
    """
    mkdir -p star_index

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SA
    touch star_index/SAindex
    touch star_index/Genome
    touch star_index/chrName.txt
    touch star_index/chrLength.txt
    touch star_index/chrStart.txt
    touch star_index/chrNameLength.txt
    touch star_index/genomeParameters.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11b
    END_VERSIONS
    """
}
