process TAMA_MERGE {
    tag "${cohort_id}"
    label 'process_high_memory'
    conda 'bioconda::gs-tama=1.0.3'
    container "${params.tama_container ?: (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0' :
        'quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0')}"

    input:
    // Shard validation deliberately emits the same basename for every BED.
    // Stage each collection member below a numbered directory so a per-accession
    // merge can accept all shard outputs without Nextflow filename collisions.
    tuple val(cohort_id), path(beds, stageAs: 'bed??/*')

    output:
    path '*_merged.bed', emit: bed
    path '*_gene_report.txt', emit: gene_report
    path '*_merge.txt', emit: merge_report
    path '*_trans_report.txt', emit: trans_report
    path 'merge_filelist.tsv', emit: filelist
    path 'versions.yml', emit: versions

    script:
    """
    prepare_tama_merge_filelist.py merge_filelist.tsv ${beds}
    tama_merge.py -f merge_filelist.tsv -p ${cohort_id}.merged -e ${params.tama_end_mode} -d merge_dup
    mv ${cohort_id}.merged.bed ${cohort_id}_merged.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tama_merge: \$(tama_merge.py -v 2>&1 | tail -n1)
    END_VERSIONS
    """

    stub:
    """
    printf 'stub\\t0\\t4\\tmerged.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${cohort_id}_merged.bed
    touch ${cohort_id}.merged_gene_report.txt ${cohort_id}.merged_merge.txt ${cohort_id}.merged_trans_report.txt
    printf 'stub\\tstub\\n' > merge_filelist.tsv
    printf '"%s":\\n    tama_merge: stub\\n' '${task.process}' > versions.yml
    """
}
