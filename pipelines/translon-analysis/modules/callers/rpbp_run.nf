process RUN_RPBP {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0'
    input:
    tuple val(meta), path(config)
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_rpbp ?: ''
    """
    cp -r \$(dirname ${config}) raw
    run-all-rpbp-instances ${config} --profiles-only --num-cpus ${task.cpus ?: 1} --logging-level INFO ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'Rp-Bp produced no output' >&2; exit 1; }
    printf '"%s":\n    Rp-Bp: 3.0.1\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf 'chr1\t300\t420\tORF_RP1\t+\t8.5\n' > raw/rpbp.bed
    printf '"stub":\n    Rp-Bp: stub\n' > versions.yml
    """
}
