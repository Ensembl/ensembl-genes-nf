// GENBLAST — protein-to-genome alignment using GenBlastG
// Custom module: no nf-core equivalent available.
// Mirrors HiveGenBlast invocation from IGTR_subpipeline.pm.

process GENBLAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::genblast=1.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genblast:1.0.4--h9f5acd7_3' :
        'biocontainers/genblast:1.0.4--h9f5acd7_3' }"

    input:
    tuple val(meta), path(proteins)
    path  genome                      // softmasked genome FASTA

    output:
    tuple val(meta), path("*.gff"),   emit: gff
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: '-P blast -gff -e 1 -c 0.8 -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000'
    """
    genblast \\
        -p genblastg \\
        -q ${proteins} \\
        -t ${genome} \\
        -o ${prefix} \\
        -g T \\
        -pid \\
        -r ${params.genblast_max_rank} \\
        -num_threads ${task.cpus} \\
        ${args}

    # GenBlast writes <prefix>.gff; rename to include sample id
    [ -f ${prefix}.gff ] && mv ${prefix}.gff ${prefix}.genblast.gff || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genblast: \$(genblast 2>&1 | grep -i version | head -1 | sed 's/.*version //' || echo '1.0.4')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    cat > ${prefix}.genblast.gff << 'STUB_GFF'
    ##gff-version 3
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=STUB-R1-1-A1;Name=STUB;PID=85.50;Coverage=98.30;Note=stub
    chr1\tgenBlastG\tcoding_exon\t1000\t2000\t.\t+\t0\tID=STUB-R1-1-A1-E1;Parent=STUB-R1-1-A1
    STUB_GFF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genblast: 1.0.4
    END_VERSIONS
    """
}
