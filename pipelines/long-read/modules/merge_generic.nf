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
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done > merged.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merged.gtf; then state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope '${scope}' --status "\${state}"
    printf 'backend\\tscope\\tstatus\\n%s\\t%s\\t%s\\n' '${backend}' '${scope}' "\${state}" > merge_status.tsv
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
