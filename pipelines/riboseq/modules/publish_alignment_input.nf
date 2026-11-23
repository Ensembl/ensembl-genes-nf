/* Preserve the exact collapsed FASTA supplied to STAR alongside its BAMs.
 * FASTA headers carry the collapsed-sequence multiplicities needed to
 * reconstruct read-weighted alignment fates from STAR's NH-tagged BAM. */

process PUBLISH_ALIGNMENT_INPUT {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/star_align/input", mode: 'copy', saveAs: {
        "${meta.id}.collapsed.fa.gz"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('published.collapsed.fa.gz'), emit: collapsed_fasta

    script:
    def read_command = reads.name.endsWith('.gz') ? 'zcat' : 'cat'
    """
    $read_command "${reads}" | gzip -c > published.collapsed.fa.gz
    """

    stub:
    """
    printf '>stub_xCOUNT1\\nACGT\\n' | gzip -c > published.collapsed.fa.gz
    """
}
