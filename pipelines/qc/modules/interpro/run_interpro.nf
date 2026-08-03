process INTERPRO_RUN {

    tag { meta.id }
    publishDir "${params.outdir}/qc/interpro", mode: 'copy', overwrite: true

    container 'docker://quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1'

    containerOptions {
        data_file_path ?
            "-B ${data_file_path.resolve()}:/opt/interproscan/data" :
             ""}
    input:
        tuple val(meta), path(protein)
        val database
        val data_file_path


    output:
        tuple val(meta), path("${meta}.tsv"), emit: stats_txt

    script:
        """
        mkdir -p tmp
        interproscan.sh \\
              -i ${protein} \\
               -b ${meta} \\
               -appl ${database}
              -T tmp \\
              -f tsv
              -cpu ${task.cpus}
        """
}
