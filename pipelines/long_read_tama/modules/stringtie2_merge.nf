process STRINGTIE2_MERGE {
    tag "${backend}:${accession}:stringtie2-merge"
    label 'process_high_memory'
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0'

    input:
    tuple val(meta), val(backend), val(accession), path(gtfs, stageAs: 'gtf??/*')

    output:
    tuple val(meta), val(backend), val(accession), path("${backend}_${accession}_merged.gtf"), emit: gtf
    path 'merge_filelist.tsv', emit: filelist
    path 'merge_filelist.sha256', emit: filelist_checksum
    path 'versions.yml', emit: versions

    script:
    def args = params.stringtie2_merge_args ?: ''
    """
    printf '%s\\n' ${gtfs} | LC_ALL=C sort > merge_filelist.tsv
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    stringtie --merge -p ${task.cpus} ${args} -o ${backend}_${accession}_merged.gtf \\
        \$(cat merge_filelist.tsv)
    test -s ${backend}_${accession}_merged.gtf || { echo 'StringTie2 merge produced an empty GTF' >&2; exit 1; }
    printf '"%s":\\n    stringtie: 2.2.3\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\tStringTie\\texon\\t1\\t4\\t.\\t+\\t.\\tgene_id "${accession}.g1"; transcript_id "${accession}.t1";\\n' > ${backend}_${accession}_merged.gtf
    printf 'stub.gtf\\n' > merge_filelist.tsv
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    printf '"%s":\\n    stringtie: 2.2.3-stub\\n' '${task.process}' > versions.yml
    """
}
