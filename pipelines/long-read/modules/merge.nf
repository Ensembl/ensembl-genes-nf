process TAMA_MERGE_ACCESSION {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_run/${meta.run_id}/models/${backend}/merge/${filename}" }
    tag "${meta.id}:tama:run"
    input:
    tuple val(meta), val(backend), path(gtfs)
    output:
    tuple val(meta), val(backend), val('run'), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    """
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done > merged.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merged.gtf; then state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope run --status "\${state}"
    printf 'backend\\tscope\\tstatus\\n%s\\trun\\t%s\\n' '${backend}' "\${state}" > merge_status.tsv
    printf '"%s":\\n    merge: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf; printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"run","stub":true}\\n' '${backend}' > merge_manifest.json; printf 'backend\\tscope\\tstatus\\n%s\\trun\\tCOMPLETE\\n' '${backend}' > merge_status.tsv; printf '"%s":\\n    merge: stub\\n' '${task.process}' > versions.yml
    """
}

process TAMA_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_cohort/${meta.cohort_id}/models/${backend}/merge/${filename}" }
    tag "${meta.id}:tama:cohort"
    input:
    tuple val(meta), val(backend), path(gtfs)
    output:
    tuple val(meta), val(backend), val('cohort'), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    """
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done > merged.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merged.gtf; then state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope cohort --status "\${state}"; printf 'backend\\tscope\\tstatus\\n%s\\tcohort\\t%s\\n' '${backend}' "\${state}" > merge_status.tsv; printf '"%s":\\n    merge: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf; printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"cohort","stub":true}\\n' '${backend}' > merge_manifest.json; printf 'backend\\tscope\\tstatus\\n%s\\tcohort\\tCOMPLETE\\n' '${backend}' > merge_status.tsv; printf '"%s":\\n    merge: stub\\n' '${task.process}' > versions.yml
    """
}

process STRINGTIE2_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_${scope == 'run' ? 'run/' + meta.run_id : 'cohort/' + meta.cohort_id}/models/${backend}/merge/${filename}" }
    // audit container 'quay.io/biocontainers/stringtie@sha256:84207eb013dbcff90928d8ddc75902580d3885ce8edaa36b8fab7665a3dc2e3c'
    container params.use_containers ? params.tool_containers.stringtie2 : null
    tag "${meta.id}:stringtie2:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtfs)
    val merge_args
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done > merge_input.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merge_input.gtf; then stringtie --merge ${args} -o merged.gtf ${gtfs}; state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope '${scope}' --status "\${state}"; printf 'backend\\tscope\\tstatus\\n%s\\t%s\\t%s\\n' '${backend}' '${scope}' "\${state}" > merge_status.tsv; printf '"%s":\\n    stringtie: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf; printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"%s","stub":true}\\n' '${backend}' '${scope}' > merge_manifest.json; printf 'backend\\tscope\\tstatus\\n%s\\t%s\\tCOMPLETE\\n' '${backend}' '${scope}' > merge_status.tsv; printf '"%s":\\n    stringtie: stub\\n' '${task.process}' > versions.yml
    """
}

process STRINGTIE3_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_${scope == 'run' ? 'run/' + meta.run_id : 'cohort/' + meta.cohort_id}/models/${backend}/merge/${filename}" }
    // audit container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34d25ca28cb224a73658053a1f31'
    container params.use_containers ? params.tool_containers.stringtie3 : null
    tag "${meta.id}:stringtie3:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(gtfs)
    val merge_args
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done > merge_input.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merge_input.gtf; then stringtie --merge ${args} -o merged.gtf ${gtfs}; state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope '${scope}' --status "\${state}"; printf 'backend\\tscope\\tstatus\\n%s\\t%s\\t%s\\n' '${backend}' '${scope}' "\${state}" > merge_status.tsv; printf '"%s":\\n    stringtie: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf; printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"%s","stub":true}\\n' '${backend}' '${scope}' > merge_manifest.json; printf 'backend\\tscope\\tstatus\\n%s\\t%s\\tCOMPLETE\\n' '${backend}' '${scope}' > merge_status.tsv; printf '"%s":\\n    stringtie: stub\\n' '${task.process}' > versions.yml
    """
}

process TMERGE_MERGE_ACCESSION {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_run/${meta.run_id}/models/${backend}/merge/${filename}" }
    tag "${meta.id}:tmerge:run"
    input:
    tuple val(meta), val(backend), path(gtfs)
    val merge_args
    output:
    tuple val(meta), val(backend), val('run'), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done | LC_ALL=C sort -k1,1 -k4,4n -k5,5n > merge_input.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merge_input.gtf; then perl ${projectDir}/bin/tmerge ${args} merge_input.gtf > merged.gtf; state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope run --status "\${state}"; printf 'backend\\tscope\\tstatus\\n%s\\trun\\t%s\\n' '${backend}' "\${state}" > merge_status.tsv; printf '"%s":\\n    tmerge: vendored-v1.0-8b4d6e7c1c94955931946081476e326c4cece161\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf; printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"run","stub":true}\\n' '${backend}' > merge_manifest.json; printf 'backend\\tscope\\tstatus\\n%s\\trun\\tCOMPLETE\\n' '${backend}' > merge_status.tsv; printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}

process TMERGE_NATIVE_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merge_*', saveAs: { filename -> "per_cohort/${meta.cohort_id}/models/${backend}/merge/${filename}" }
    tag "${meta.id}:tmerge:cohort"
    input:
    tuple val(meta), val(backend), path(gtfs)
    val merge_args
    output:
    tuple val(meta), val(backend), val('cohort'), path('merged.gtf'), path('merge_filelist.txt'), path('merge_filelist.sha256'), emit: merged
    path 'merge_manifest.json', emit: manifest
    path 'merge_status.tsv', emit: status
    path 'merge_skipped_inputs.tsv', emit: skipped
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: merge_args
    """
    for input_gtf in \$(printf '%s\\n' ${gtfs} | sort); do cat "\${input_gtf}"; done | LC_ALL=C sort -k1,1 -k4,4n -k5,5n > merge_input.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' merge_input.gtf; then perl ${projectDir}/bin/tmerge ${args} merge_input.gtf > merged.gtf; state=COMPLETE; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv; else printf '# EMPTY: no native exon records\\n' > merged.gtf; state=EMPTY; printf 'input\\tstatus\\n%s\\tSKIPPED_EMPTY\\n' '${backend}' > merge_skipped_inputs.tsv; fi
    python3 ${projectDir}/bin/merge_filelist.py --inputs ${gtfs} --file-list merge_filelist.txt --checksum merge_filelist.sha256 --manifest merge_manifest.json --backend '${backend}' --scope cohort --status "\${state}"; printf 'backend\\tscope\\tstatus\\n%s\\tcohort\\t%s\\n' '${backend}' "\${state}" > merge_status.tsv; printf '"%s":\\n    tmerge: vendored-v1.0-8b4d6e7c1c94955931946081476e326c4cece161\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cat ${gtfs} > merged.gtf; printf 'stub\\n' > merge_filelist.txt; printf 'stub\\n' > merge_filelist.sha256; printf 'input\\tstatus\\n' > merge_skipped_inputs.tsv
    printf '{"backend":"%s","scope":"cohort","stub":true}\\n' '${backend}' > merge_manifest.json; printf 'backend\\tscope\\tstatus\\n%s\\tcohort\\tCOMPLETE\\n' '${backend}' > merge_status.tsv; printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
