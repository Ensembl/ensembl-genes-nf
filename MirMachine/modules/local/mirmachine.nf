
process MIRMACHINE {
    label 'mirMachine'

    tag "${meta.id}"
    publishDir "${params.outdir}/mirmachine/${meta.id}", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    maxForks 10

    input:
    tuple val(meta), path(fasta), val(node), val(model)

    output:
    tuple val(meta), path("results/predictions/*"), emit: predictions
    tuple val(meta), path("${meta.id}_*.heatmap.csv"), emit: csv
    path "${meta.id}_mirmachine.log", emit: log

    script:
    def species = meta.species.replace(" ", "_")
    """
    MirMachine.py --node ${node} \
                  --species ${species} \
                  --genome ${fasta} \
                  --model ${model} \
                  --cpu ${task.cpus}  2> ${meta.id}_mirmachine.log

    cp results/predictions/heatmap/*.heatmap.csv  ./${meta.id}_${species}.heatmap.csv
    """
}