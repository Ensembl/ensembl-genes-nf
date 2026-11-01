process MERGE_BACKEND_MODELS {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_${scope == 'run' ? 'run/' + meta.run_id : 'cohort/' + meta.cohort_id}/models/${backend}/merge/${filename}" }
    tag "${meta.id}:${backend}"
    input:
    tuple val(meta), val(backend), path(gtfs), val(scope)
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${gtfs} --output merged.gtf --backend '${backend}' --status merge_gtf.status.tsv --skipped merge_skipped_inputs.tsv --merge-status merge_status.tsv --scope '${scope}'
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope '${scope}' --status-file merge_gtf.status.tsv
    printf '"%s":\\n    merge: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf
    printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"%s","stub":true}\\n' '${backend}' '${scope}' > merge_manifest.json
    printf 'backend\\tscope\\tstatus\\n%s\\t%s\\tCOMPLETE\\n' '${backend}' '${scope}' > merge_status.tsv
    printf '"%s":\\n    merge: stub\\n' '${task.process}' > versions.yml
    """
}
