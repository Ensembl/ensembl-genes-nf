process INTERPRO_RUN {

    tag { meta.id }
    publishDir "${params.outdir}/qc/interpro", mode: 'copy', overwrite: true

    // InterProScan 5 image; Singularity pulls this Docker image under -profile slurm.
    container 'docker://interpro/interproscan:5.78-109.0'

    containerOptions {
        data_file_path ? "-B ${file(data_file_path).resolve()}:/opt/interproscan/data:ro" : ""
    }

    input:
        tuple val(meta), path(protein)
        val database
        val data_file_path

    output:
        tuple val(meta), path("${meta.sample ?: meta.id ?: protein.simpleName}.tsv"), emit: stats_txt
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def stem = meta.sample ?: meta.id ?: protein.simpleName
        def args = task.ext.args ?: ''
        def interproscan = '/opt/interproscan/interproscan.sh'
        """
        mkdir -p tmp
        ${interproscan} \\
              -i ${protein} \\
               -b ${stem} \\
               -appl ${database} \\
              -T tmp \\
              -f tsv \\
              -cpu ${task.cpus} \\
              ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            interproscan: \$(${interproscan} -version 2>&1 | head -n 1 | sed 's/^.*InterProScan //')
        END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: protein.simpleName
        """
        touch ${stem}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            interproscan: 5.78-109.0
        END_VERSIONS
        """
}
