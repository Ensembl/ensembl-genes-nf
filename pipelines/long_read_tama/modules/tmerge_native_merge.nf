process TMERGE_NATIVE_MERGE {
    tag "${backend}:${accession}:tmerge-merge"
    label 'process_high_memory'
    container 'community.wave.seqera.io/library/pip_pyfaidx_setuptools_six_pruned:b4fcc6bbecfa2bea'

    input:
    tuple val(meta), val(backend), val(accession), path(gtfs, stageAs: 'gtf??/*')

    output:
    tuple val(meta), val(backend), val(accession), path("${backend}_${accession}_merged.gtf"), emit: gtf
    path 'merge_filelist.tsv', emit: filelist
    path 'merge_filelist.sha256', emit: filelist_checksum
    path 'versions.yml', emit: versions

    script:
    """
    printf '%s\\n' ${gtfs} | LC_ALL=C sort > merge_filelist.tsv
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    cat \$(cat merge_filelist.tsv) > tmerge_native_input.gtf
    test -s tmerge_native_input.gtf || { echo 'tmerge native merge received no GTF records' >&2; exit 1; }
    tmerge ${params.tmerge_args ?: ''} --input tmerge_native_input.gtf --output ${backend}_${accession}_merged.gtf
    test -s ${backend}_${accession}_merged.gtf || { echo 'tmerge native merge produced an empty GTF' >&2; exit 1; }
    printf '"%s":\\n    tmerge: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\ttmerge\\texon\\t1\\t4\\t.\\t+\\t.\\ttranscript_id "${accession}.t1";\\n' > ${backend}_${accession}_merged.gtf
    printf 'stub.gtf\\n' > merge_filelist.tsv
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
