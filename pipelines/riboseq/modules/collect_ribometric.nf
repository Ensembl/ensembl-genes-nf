process COLLECT_RIBOMETRIC {
    tag "${meta.id}"
    label "process_medium"

    container "python:3.10-slim"

    input:
    val run_id
    tuple val(meta), path(json), path(csv), path(offsets)

    output:
    path "ribometric_metrics.done", emit: done

    script:
    def db = params.metrics_db ?: "${params.outdir}/pipeline_info/qc.duckdb"
    """
    python3 -m pip install -q --no-cache-dir duckdb pandas >/dev/null 2>&1 || true
    python3 $projectDir/bin/collect_ribometric.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --json ${json} \
      --csv ${csv} \
      --offsets ${offsets}
    touch ribometric_metrics.done
    """
}
