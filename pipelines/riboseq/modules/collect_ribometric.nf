process COLLECT_RIBOMETRIC {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.10 conda-forge::duckdb conda-forge::pandas"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    input:
    val run_id
    tuple val(meta), path(json), path(csv), path(offsets)

    output:
    tuple val(meta), path("ribometric_metrics.done"), emit: done

    script:
    def db = params.metrics_db ?: "${params.outdir}/pipeline_info/qc.duckdb"
    """
    collect_ribometric.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --json ${json} \
      --csv ${csv} \
      --offsets ${offsets}
    touch ribometric_metrics.done
    """
}
