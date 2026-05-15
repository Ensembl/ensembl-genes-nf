
process MIRMACHINE {
    label 'mirMachine'

    tag "${meta.id}"
    publishDir "${params.outdir}/mirmachine/${meta.id}", mode: 'copy'

    // errorStrategy { task.attempt <= 0 ? 'retry' : 'ignore' }
    maxRetries 3

    memory { 20.GB * task.attempt }
    
    time 24.h
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    
    maxRetries 5

    maxForks 100

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/predictions/*"), emit: predictions
    tuple val(meta), path("${meta.id}_*.heatmap.csv"), emit: csv
    path "${meta.id}_mirmachine.log", emit: log

    script:
    def species = meta.species.replace(" ", "_")
    def long_hairpin = params.long_hairpin ? "--long" : ""
    """
    # Set Snakemake cache directory to current working directory
    export HOME=\$PWD

    MirMachine.py --species ${species} \
                  --genome ${fasta} \
                  --family ${params.family} \
                  --model ${params.model} \
                  --evalue ${params.evalue} \
                  --cpu ${task.cpus} \
                  ${long_hairpin} 2> ${meta.id}_mirmachine.log

    cp results/predictions/heatmap/*.heatmap.csv  ./${meta.id}_${species}.heatmap.csv
    """
}
