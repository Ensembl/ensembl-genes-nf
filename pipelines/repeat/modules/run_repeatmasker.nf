#!/usr/bin/env nextflow
/*
This process runs RepeatMasker on a given genome file using a 
specified RepeatModeler library. It uses the RepeatMasker tool
to perform the analysis and generates GTF files containing the identified 
repetitive regions. The output GTF files are saved in the "repeatmasker" 
directory under the output directory for the given GCA accession. 
The process also generates a versions.yml file containing the version 
of RepeatMasker used.
*/
process RUN_REPEATMASKER {
    label "python"
    tag "${meta.gca}:genome"

    publishDir "${params.outdir}/${meta.gca}/repeatmasker/", pattern: "repeatmasker_output/*.fa", mode: "copy"

    input:
    tuple val(meta), val(library_file)

    output:
    tuple val(meta), path("*genome.masked.fa"), emit: repeatmasker_softmasked_out
    tuple val(meta), path("*.gtf"), emit: repeatmasker_out
    path "versions.yml", emit: versions_file

    script:
    def library = meta.repeatmasker_library?.trim() ? meta.repeatmasker_library : library_file
    LD_LIBRARY_PATH='/hps/software/users/ensembl/genebuild/shared/libnsl/lib:/hps/software/users/ensembl/genebuild/shared/libnsl/libtirpc'
    // /opt/linuxbrew/bin/RepeatMasker -help
    """
    export LD_LIBRARY_PATH=/hps/software/users/ensembl/genebuild/shared/libnsl/lib:/hps/software/users/ensembl/genebuild/shared/libnsl/libtirpc:$LD_LIBRARY_PATH
    ls /hps/software/users/ensembl/genebuild/shared/libnsl/lib/libnsl.so.2
    run_repeatmasker  --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --repeatmasker_bin /opt/linuxbrew/bin/RepeatMasker \
                    --library ${library} \
                    --repeatmasker_engine ${params.engine_repeatmasker} \
                    --num_threads ${task.cpus} \
                    --bedtools_bin /opt/linuxbrew/bin/bedtools
    mv repeatmasker_output/annotation.gtf ${meta.gca}_repeatmasker.gtf 
    mv repeatmasker_output/genome.masked.fa ${meta.gca}_genome.masked.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -version 2>&1 | head -n 1 | sed 's/.*RepeatMasker version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}
