process TAMA_MERGE_ACCESSION {
    tag "${backend}:${accession}"
    label 'process_high_memory'
    conda 'bioconda::gs-tama=1.0.3'
    container "${params.tama_container ?: (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0' :
        'quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0')}"

    input:
    tuple val(meta), val(backend), val(accession), path(beds, stageAs: 'bed??/*')

    output:
    tuple val(meta), val(backend), val(accession), path("${backend}_${accession}_merged.bed"), emit: bed
    path 'merge_filelist.tsv', emit: filelist
    path 'merge_filelist.sha256', emit: filelist_checksum
    path 'accession_merge_report.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    prepare_tama_merge_filelist.py merge_filelist.tsv ${beds} \\
        --seq-type ${params.tama_cap_mode} --priority-rank 1,1,1
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    tama_merge.py -f merge_filelist.tsv -p ${accession}.merged -e ${params.tama_end_mode} -d merge_dup
    test -s ${accession}.merged.bed || { echo 'TAMA accession merge produced an empty BED' >&2; exit 1; }
    mv ${accession}.merged.bed ${backend}_${accession}_merged.bed
    printf 'backend\\taccession\\tstatus\\tinputs\\n${backend}\\t${accession}\\tSUCCESS\\t%s\\n' "$(find . -path './bed??/*' -name '*.bed' | wc -l | tr -d ' ')" > accession_merge_report.tsv
    printf '"%s":\\n    tama_merge: \$(tama_merge.py -v 2>&1 | tail -n1)\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\t0\\t4\\t${backend}.${accession}.merged.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${backend}_${accession}_merged.bed
    printf 'stub.bed\\tno_cap\\t1,1,1\\t${backend}\\n' > merge_filelist.tsv
    sha256sum merge_filelist.tsv > merge_filelist.sha256
    printf 'backend\\taccession\\tstatus\\tinputs\\n${backend}\\t${accession}\\tSUCCESS\\t1\\n' > accession_merge_report.tsv
    printf '"%s":\\n    tama_merge: stub\\n' '${task.process}' > versions.yml
    """
}
