
process CMSEARCH {

    publishDir "${params.outdir}/cmsearch_output", mode: copy

    input:
        

    output:
        

    script:
        """
        cmsearch -h
        """
}
