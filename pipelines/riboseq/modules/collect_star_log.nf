process COLLECT_STAR_LOG {
    tag "${meta.id}"
    label "process_medium"

    container "python:3.10"

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
    python3 -m pip install -q --no-cache-dir duckdb pandas >/dev/null 2>&1 || true
    python3 $projectDir/bin/collect_star_log.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --study-id ${meta.study_id ?: 'unknown'} \
      --log ${star_log}
    touch star_metrics.done
    """
}
