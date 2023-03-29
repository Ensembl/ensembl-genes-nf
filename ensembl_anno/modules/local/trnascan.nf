
process TRNASCAN {

    publishDir "${params.outdir}/trnascan_output", mode: copy

    input:
        

    output:
        

    script:
        // """
        // tRNAscan-SE -o -f -H -q --detail -Q
        """
        tRNAscan-SE -h
        """
}
