process COLLECT_GETRPF_CLEAN {
    tag "${meta.id}"

    container "python:3.10"

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
    python3 -m pip install -q --no-cache-dir duckdb >/dev/null 2>&1 || true
    if [ ! -s ${report} ] || [ ! -s ${checks} ]; then
      touch getrpf_metrics.done
      exit 0
    fi
    python3 $projectDir/bin/collect_getrpf_clean.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --report ${report} \
      --checks ${checks}
    touch getrpf_metrics.done
    """
}
