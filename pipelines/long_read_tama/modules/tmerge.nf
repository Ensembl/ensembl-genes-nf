process TMERGE {
    tag "${cohort_id}"
    label 'process_high_memory'
    // This Python tmerge image uses explicit --input/--output arguments.
    container 'community.wave.seqera.io/library/pip_pyfaidx_setuptools_six_pruned:b4fcc6bbecfa2bea'

    input:
    tuple val(cohort_id), path(beds)

    output:
    path '*_merged.bed', emit: bed
    path '*_gene_report.txt', emit: gene_report
    path '*_merge.txt', emit: merge_report
    path '*_trans_report.txt', emit: trans_report
    path 'merge_filelist.tsv', emit: filelist
    path 'versions.yml', emit: versions

    script:
    """
    # tmerge requires exon-only, coordinate-sorted GTF input with a unique
    # transcript_id for every input model.
    bed12_to_gtf.py tmerge_input.gtf ${beds}
    tmerge ${params.tmerge_args} --input tmerge_input.gtf --output tmerge_output.gtf
    test -s tmerge_output.gtf || { echo 'tmerge produced an empty GTF' >&2; exit 1; }
    gtf_to_bed12.py tmerge_output.gtf ${cohort_id}_merged.bed ${cohort_id}
    test -s ${cohort_id}_merged.bed || { echo 'tmerge BED conversion produced no models' >&2; exit 1; }
    printf 'tmerge_input.gtf\\n' > merge_filelist.tsv
    touch ${cohort_id}_merged_gene_report.txt ${cohort_id}_merged_merge.txt ${cohort_id}_merged_trans_report.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tmerge: \$(tmerge --version 2>&1 | tail -n1 || true)
    END_VERSIONS
    """

    stub:
    """
    printf 'stub\\t0\\t4\\tmerged.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${cohort_id}_merged.bed
    touch ${cohort_id}_merged_gene_report.txt ${cohort_id}_merged_merge.txt ${cohort_id}_merged_trans_report.txt
    printf 'tmerge_input.gtf\\n' > merge_filelist.tsv
    printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
