// CMSEARCH — search Rfam covariance models against genome slices.
// Matches nf-core/infernal/cmsearch interface; uses the same container.
// Outputs --tblout (tabular) which is parsed by filter_ncrna.py.

process CMSEARCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321hece9b99_0' :
        'biocontainers/infernal:1.1.5--pl5321hece9b99_0' }"

    input:
    tuple val(meta), path(cmfile), path(seqdb)

    output:
    tuple val(meta), path("*.tblout"), emit: tblout
    tuple val(meta), path("*.txt"),    emit: output
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: '--rfam --nohmmonly --cut_ga'
    """
    cmsearch \\
        ${args} \\
        --cpu ${task.cpus} \\
        --tblout ${prefix}.tblout \\
        ${cmfile} \\
        ${seqdb} \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmsearch -h 2>&1 | grep INFERNAL | sed 's/.*INFERNAL //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    cat > ${prefix}.tblout << 'STUB_TBLOUT'
    # cmsearch stub output
    chr1\t-\tmir-21\tRF00001\tcm\t1\t88\t1000\t1088\t+\t-\t1\t0.46\t0.0\t87.3\t1.2e-16\t!\tmiRNA stub
    STUB_TBLOUT
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: 1.1.5
    END_VERSIONS
    """
}
