
process MIRMACHINE {
    label 'mirMachine'

    tag "${meta.id}"
    publishDir "${params.outdir}/mirmachine/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(fasta), val(node), val(model)

    output:
    tuple val(meta), path("results/predictions/*"), emit: predictions
    path "${meta.id}_mirmachine.log", emit: log

    script:
    def species = meta.species.replace(" ", "_")
    """
    MirMachine.py --node ${node} \
                  --species ${species} \
                  --genome ${fasta} \
                  --model ${model} \
                  --cpu ${task.cpus} 2> ${meta.id}_mirmachine.log
    """
}