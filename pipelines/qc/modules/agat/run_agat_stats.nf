process AGAT_RUN_STATS {

    tag { meta.id }
    publishDir "${params.outdir}/qc/agat", mode: 'copy', overwrite: true,
        pattern: "*_agat_stats.txt"
    container 'docker://quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'

    // AGAT reads feature_levels.yaml from this path. The pipeline owns the
    // configuration so users do not need to provide a host path.
    containerOptions {
        "-B ${file(feature_levels_yaml)}:/usr/local/lib/perl5/site_perl/auto/share/dist/AGAT/feature_levels.yaml:ro"
    }

    input:
        tuple val(meta), path(gff3)
        val feature_levels_yaml

    output:
        tuple val(meta), path("*_agat_stats.txt"), emit: stats_txt
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def stem = meta.sample ?: meta.id ?: gff3.simpleName
        def args = task.ext.args ?: ''

        """
            agat_sp_statistics.pl \\
              --gff ${gff3} \\
              -o ${stem}_agat_stats.txt \\
              --cpu 0 \\
              --verbose 3 \\
              ${args}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                agat: \$(agat_sp_statistics.pl --version 2>&1 | head -n 1 | sed 's/^.*AGAT //')
            END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: gff3.simpleName
        """
        touch ${stem}_agat_stats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            agat: 1.7.0
        END_VERSIONS
        """
}
