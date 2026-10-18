process PREP_DIAMOND_DB {
    tag { 'reference' }
    label 'process_medium'
    container 'community.wave.seqera.io/library/diamond:2.1.24--61a5af76160d103f'
    publishDir "${params.outdir}/qc/diamond", mode: 'copy', overwrite: true,
        pattern: 'reference.dmnd'

    input:
        path ref_protein_faa

    output:
        path 'reference.dmnd', emit: db
        path 'versions.yml', emit: versions

    script:
        """
        diamond makedb --in ${ref_protein_faa} --db reference

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: \$(diamond version | sed 's/^diamond version //')
        END_VERSIONS
        """

    stub:
        """
        touch reference.dmnd

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: 2.1.24
        END_VERSIONS
        """
}
