
process LIST_MIRMACHINE_CLADES {

    label 'mirMachine'
    
    output:
        path "*txt"

    script:
        """
        MirMachine.py --print-all-nodes > clades.txt
        """
    
}


process FORMAT_CLADES {

    input:
        file clades

    output:
        file "formatted_clades.txt"

    shell:
        """ 
        awk 'NR>1 {for (i=1; i<=NF; i++) print \$i}' ${clades} > formatted_clades.txt
        """

}

process MATCH_CLADE {

    //label 'mirMachine'
    errorStrategy 'retry'
    maxRetries 3

    input:
        val(species)
        val(accession)
        file clade_file

    output:
        stdout

    script:
        """
        python $projectDir/scripts/match_clade.py -s "${species}" -c $clade_file --output stdout
        """

}
