// AUGUSTUS
// Run Augustus ab initio gene prediction on a genome chunk.
//
// Species model: built-in name (e.g. 'human', 'fly') or path to a custom
// config directory (set via params.augustus_config_path).

process AUGUSTUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::augustus=3.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.5.0--pl5321h700735d_3' :
        'biocontainers/augustus:3.5.0--pl5321h700735d_3' }"

    input:
    tuple val(meta), path(genome_chunk)
    val   species

    output:
    tuple val(meta), path("*.aug.gff"), emit: gff
    path  "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: meta.id
    def args       = task.ext.args   ?: ''
    def config_opt = params.augustus_config_path
        ? "--AUGUSTUS_CONFIG_PATH=${params.augustus_config_path}"
        : ''
    def hints_opt  = meta.hints_bam ? "--hintsfile=${meta.hints_bam} --extrinsicCfgFile=${params.extrinsic_cfg}" : ''
    """
    augustus \\
        --species=${species} \\
        --gff3=on \\
        --softmasking=1 \\
        --outfile=${prefix}.aug.gff \\
        ${config_opt} \\
        ${hints_opt} \\
        ${args} \\
        ${genome_chunk}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augustus: \$(augustus --version 2>&1 | head -1 | grep -oP '[\\d.]+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\\n' > ${prefix}.aug.gff
    printf 'chr1\\taugustus\\tgene\\t1000\\t5000\\t.\\t+\\t.\\tID=g1\\n' >> ${prefix}.aug.gff
    printf 'chr1\\taugustus\\tmRNA\\t1000\\t5000\\t.\\t+\\t.\\tID=g1.t1;Parent=g1\\n' >> ${prefix}.aug.gff
    printf 'chr1\\taugustus\\texon\\t1000\\t2000\\t.\\t+\\t.\\tID=g1.t1.exon1;Parent=g1.t1\\n' >> ${prefix}.aug.gff
    printf 'chr1\\taugustus\\texon\\t3000\\t5000\\t.\\t+\\t.\\tID=g1.t1.exon2;Parent=g1.t1\\n' >> ${prefix}.aug.gff
    printf 'chr1\\taugustus\\tCDS\\t1000\\t2000\\t.\\t+\\t0\\tID=g1.t1.cds1;Parent=g1.t1\\n' >> ${prefix}.aug.gff
    printf 'chr1\\taugustus\\tCDS\\t3000\\t5000\\t.\\t+\\t0\\tID=g1.t1.cds2;Parent=g1.t1\\n' >> ${prefix}.aug.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augustus: 3.5.0
    END_VERSIONS
    """
}
