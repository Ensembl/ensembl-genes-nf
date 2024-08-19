
process MIRMACHINE {
    //container "${params.container}"
    def formatted_species_name = params.species.replace(" ","_")
    publishDir "${params.outdir}/mirmachine/${formatted_species_name}_${params.accession}", mode: 'copy'

    //label 'mirMachine'

    input:
        val(species)
        val(node)
        val(model)
        file fasta

    output:
        path "results/predictions/*"

    script:
        """
        singularity exec ${params.container} /bin/bash -c "
            MirMachine.py --cpu ${task.cpus} --node ${node} --model ${model} --species  ${species} --genome ${fasta} > log_file 2>> err_file
        "
        """
}
