#!/usr/bin/env nextflow
/*
This process runs the Red to identify repetitive regions in a genome file. 
It uses the Red tool to perform the analysis and generates a GTF file 
containing the identified repetitive regions. The output GTF file is 
saved in the "red" directory under the output directory for the given 
GCA accession. The process also generates a versions.yml 
file containing the version of Red used.
*/
process RUN_RED {
    label "python"
    tag "${meta.gca}:genome"
    publishDir "${params.outdir}/${meta.gca}/red/", pattern: "**/*.gtf", mode: "copy"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.gtf"), emit: red_out
    path "versions.yml", emit: versions_file    

    script:
    """
    run_red --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --red_bin ${params.red_path}
    mv red_output/annotation.gtf ${meta.gca}_red.gtf  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        red: \$(Red --version 2>&1 | head -n 1 | sed 's/.*Red version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}
