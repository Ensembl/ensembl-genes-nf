// STAR_INDEX
// Build a STAR genome index from an unmasked genome FASTA.
// The index directory is published to outdir/star_index/.

process STAR_INDEX {
    label 'process_high'

    conda "bioconda::star=2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    publishDir "${params.outdir}/star_index", mode: 'copy'

    input:
    path genome_fasta
    path gtf           // optional annotation for junction database; pass [] to skip

    output:
    path "star_index/", emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def gtf_arg  = gtf ? "--sjdbGTFfile ${gtf} --sjdbOverhang ${params.sjdb_overhang}" : ''
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        ${gtf_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/Genome star_index/SA star_index/SAindex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11b
    END_VERSIONS
    """
}
