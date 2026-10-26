process STRINGTIE3_MERGE {
    publishDir path: { accession == params.cohort_id ? "${params.outdir}/per_cohort/${params.cohort_id}/tools/stringtie3" : "${params.outdir}/per_run/${accession}/tools/stringtie3" }, mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${backend}:${accession}:stringtie3-merge"
    label 'process_high_memory'
    container 'https://depot.galaxyproject.org/singularity/stringtie:3.0.3--h29c0135_0'

    input:
    tuple val(meta), val(backend), val(accession), path(gtfs, stageAs: 'gtf??/*')

    output:
    tuple val(meta), val(backend), val(accession), path("${backend}_${accession}_merged.gtf"), emit: gtf
    path 'merge_filelist.tsv', emit: filelist
    path 'merge_skipped_empty.tsv', emit: skipped_empty
    path 'merge_filelist.sha256', emit: filelist_checksum
    path 'versions.yml', emit: versions

    script:
    def args = params.stringtie3_merge_args ?: ''
    """
    : > merge_filelist.tsv
    : > merge_skipped_empty.tsv
    for gtf in ${gtfs}; do
        if awk 'NF >= 9 && \$0 !~ /^#/ { found=1; exit } END { exit !found }' "\$gtf"; then
            printf '%s\\n' "\$gtf" >> merge_filelist.tsv
        else
            printf '%s\\n' "\$gtf" >> merge_skipped_empty.tsv
        fi
    done
    LC_ALL=C sort -o merge_filelist.tsv merge_filelist.tsv
    LC_ALL=C sort -o merge_skipped_empty.tsv merge_skipped_empty.tsv
    test -s merge_filelist.tsv || { echo 'StringTie3 merge received no GTFs with records; see merge_skipped_empty.tsv' >&2; exit 1; }
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    stringtie --merge -p ${task.cpus} ${args} -o ${backend}_${accession}_merged.gtf \\
        \$(cat merge_filelist.tsv)
    test -s ${backend}_${accession}_merged.gtf || { echo 'StringTie3 merge produced an empty GTF' >&2; exit 1; }
    printf '"%s":\\n    stringtie: 3.0.3\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\tStringTie\\texon\\t1\\t4\\t.\\t+\\t.\\tgene_id "${accession}.g1"; transcript_id "${accession}.t1";\\n' > ${backend}_${accession}_merged.gtf
    printf 'stub.gtf\\n' > merge_filelist.tsv
    : > merge_skipped_empty.tsv
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    printf '"%s":\\n    stringtie: 3.0.3-stub\\n' '${task.process}' > versions.yml
    """
}
