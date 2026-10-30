#!/usr/bin/env nextflow
/*
This process runs DustMasker to identify low-complexity regions in a genome file. 
It uses the dustmasker tool to perform the analysis and generates a GTF file 
containing the identified low-complexity regions. The output GTF file is saved 
in the "dust" directory under the output directory for the given GCA accession. 
The process also generates a versions.yml file containing the version of Dust used.
*/

process RUN_DUST {
    label "python"
    tag "${meta.gca}:genome"

    publishDir "${params.outdir}/${meta.gca}/dust/", pattern: "**/*.gtf", mode: "copy"
    input:
    val(meta)

    output:
    tuple val(meta), path("*.gtf"), emit: dust_out
    path "versions.yml", emit: versions_file

    script:
    """
    run_dust --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --dust_bin /opt/linuxbrew/bin/dustmasker \
                    --num_threads ${task.cpus}  \
                    --bedtools_bin /opt/linuxbrew/bin/bedtools
    mv dust_output/annotation.gtf ${meta.gca}_dust.gtf  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dust: \$(Dust -version 2>&1 | head -n 1 | sed 's/.*Dust version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
