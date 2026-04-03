// PREPARE_GENOME
// Decompress and index the genome FASTA with samtools faidx.
// The indexed genome is the primary output for downstream analyses.

process PREPARE_GENOME {
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path fasta_gz

    output:
    path "*.fna",     emit: fasta
    path "*.fna.fai", emit: fai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out_fasta = fasta_gz.name.replaceAll(/\.gz$/, '')
    """
    bgzip -d -c ${fasta_gz} > ${out_fasta}
    samtools faidx ${out_fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def out_fasta = fasta_gz.name.replaceAll(/\.gz$/, '')
    """
    printf '>chr1\\nATCGATCG\\n' > ${out_fasta}
    printf 'chr1\\t8\\t6\\t8\\t9\\n' > ${out_fasta}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
