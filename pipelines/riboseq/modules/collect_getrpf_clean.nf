process COLLECT_GETRPF_CLEAN {
    tag "${meta.id}"

    conda "conda-forge::python=3.10 conda-forge::duckdb"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    input:
    val run_id
    tuple val(meta), path(report)
    tuple val(meta2), path(checks)

    output:
    tuple val(meta), path("getrpf_metrics.done"), emit: done

    when:
    task.ext.when == null || task.ext.when

    script:
    def db = params.metrics_db ?: "${params.outdir}/pipeline_info/qc.duckdb"
    """
    if [ ! -s ${report} ] || [ ! -s ${checks} ]; then
      touch getrpf_metrics.done
      exit 0
    fi
    collect_getrpf_clean.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --report ${report} \
      --checks ${checks}
    touch getrpf_metrics.done
    """
}
