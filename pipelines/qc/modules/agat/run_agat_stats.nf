process AGAT_RUN_STATS {

    tag { meta.id }
    publishDir "${params.outdir}/qc/agat", mode: 'copy', overwrite: true,
        pattern: "*_agat_stats.txt"
    container 'docker://quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'

    
    containerOptions {
        feature_levels_yaml ?
            "-B ${feature_levels_yaml.resolve()}:/usr/local/lib/perl5/site_perl/auto/share/dist/AGAT/feature_levels.yaml:ro" :
             ""}
    
    input:
        tuple val(meta), path(gff3)
        val  feature_levels_yaml

    output:
        tuple val(meta), path("${meta.sample ?: meta.id ?: gff3.simpleName}_agat_stats.txt"), emit: stats_txt

    script:
        def stem     = meta.sample ?: meta.id ?: gff3.simpleName

        """
            agat_sp_statistics.pl \\
              --gff ${gff3} \\
              -o ${stem}_agat_stats.txt
        """
}
