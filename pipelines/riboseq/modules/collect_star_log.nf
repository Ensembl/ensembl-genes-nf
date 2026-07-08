process COLLECT_STAR_LOG {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.10 conda-forge::duckdb"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    input:
    val run_id
    tuple val(meta), path(star_log)

    output:
    tuple val(meta), path("star_metrics.done"), emit: done

    when:
    task.ext.when == null || task.ext.when

    script:
    def db = params.metrics_db ?: "${params.outdir}/pipeline_info/qc.duckdb"
    """
    collect_star_log.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --study-id ${meta.study_id ?: 'unknown'} \
      --log ${star_log}
    touch star_metrics.done
    """
}
