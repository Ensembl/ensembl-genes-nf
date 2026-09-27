process TMERGE_MERGE_ACCESSION {
    tag "${backend}:${accession}"
    label 'process_high_memory'
    container 'community.wave.seqera.io/library/pip_pyfaidx_setuptools_six_pruned:b4fcc6bbecfa2bea'

    input:
    tuple val(backend), val(accession), path(beds, stageAs: 'bed??/*')

    output:
    tuple val(backend), val(accession), path("${backend}_${accession}_merged.bed"), emit: bed
    path 'merge_filelist.tsv', emit: filelist
    path 'accession_merge_report.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    bed12_to_gtf.py tmerge_input.gtf ${beds}
    tmerge ${params.tmerge_args} --input tmerge_input.gtf --output tmerge_output.gtf
    test -s tmerge_output.gtf || { echo 'tmerge accession merge produced an empty GTF' >&2; exit 1; }
    gtf_to_bed12.py tmerge_output.gtf ${backend}_${accession}_merged.bed ${accession}
    test -s ${backend}_${accession}_merged.bed || { echo 'tmerge accession conversion produced no models' >&2; exit 1; }
    printf 'tmerge_input.gtf\\n' > merge_filelist.tsv
    printf 'backend\\taccession\\tstatus\\tinputs\\n${backend}\\t${accession}\\tSUCCESS\\t%s\\n' "$(find . -path './bed??/*' -name '*.bed' | wc -l | tr -d ' ')" > accession_merge_report.tsv
    printf '"%s":\\n    tmerge: \$(tmerge --version 2>&1 | tail -n1 || true)\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\t0\\t4\\t${backend}.${accession}.merged.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${backend}_${accession}_merged.bed
    printf 'tmerge_input.gtf\\n' > merge_filelist.tsv
    printf 'backend\\taccession\\tstatus\\tinputs\\n${backend}\\t${accession}\\tSUCCESS\\t1\\n' > accession_merge_report.tsv
    printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
