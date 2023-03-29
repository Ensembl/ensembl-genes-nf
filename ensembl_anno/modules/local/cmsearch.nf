
process CMSEARCH {

    publishDir "${params.outdir}/cmsearch_output", mode: copy

    input:
        

    output:
        stdout

    script:
        """
        cmsearch -h
        """
}
