
process TRNASCAN {

    // publishDir "${params.outdir}/trnascan_output", mode: copy

    // container "${projectDir}/singularity/trnascan-se%3A2.0.9--pl5321hec16e2b_3"

    input:
        

    output:
        stdout

    script:
        """
        tRNAscan-SE --help
        """
}
