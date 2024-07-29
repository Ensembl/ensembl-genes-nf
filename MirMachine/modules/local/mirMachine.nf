
process MIRMACHINE {
    //container "${params.container}"
    publishDir "${params.outdir}/mirmachine/${params.species}_${params.accession}", mode: 'copy'
    
    //label 'mirMachine'
    
    input:
        val(species)
        val(node)
        file fasta

    output:
        path "results/predictions/*"

    script:
        """
        singularity exec ${params.container} /bin/bash -c "
            MirMachine.py --cpu ${task.cpus} --node ${node} --species  ${species} --genome ${fasta} > log_file 2>> err_file
        "
        """
}
