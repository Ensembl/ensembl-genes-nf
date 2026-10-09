process ADD_METRICS_TO_REGISTRY {
    label 'python'

    tag { meta.id }
    publishDir "${params.outdir}/qc/registry", mode: 'copy', overwrite: true

    input:
        tuple val(meta), path(metrics_csv)

        tuple val(registry_host),
              val(registry_port),
              val(registry_user),
              val(registry_password),
              val(registry_db)

    output:
        tuple val(meta), path("${meta.id}.registry_metrics.done"), emit: loaded
        path 'versions.yml', emit: versions

    when:
        params.add_metrics_to_registry

    script:
        """
        metrics_to_registry.py \\
            ${metrics_csv} \\
            --assembly ${meta.id} \\
            --registry-host ${registry_host} \\
            --registry-port ${registry_port} \\
            --registry-user ${registry_user} \\
            --registry-password '${registry_password}' \\
            --registry-db ${registry_db}

        touch ${meta.id}.registry_metrics.done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | awk '{print \$2}')
            pymysql: \$(python -c 'import importlib.metadata; print(importlib.metadata.version("PyMySQL"))')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.registry_metrics.done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: stub
            pymysql: stub
        END_VERSIONS
        """
}
