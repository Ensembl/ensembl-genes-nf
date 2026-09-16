process PREDICT_LONGEST_ATG_ORFS {
    tag "${meta.id}:orfs"
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(transcripts)

    output:
    tuple val(meta), path('combined_transcripts.faa'), emit: peptides
    tuple val(meta), path('combined_orf_manifest.tsv'), emit: manifest
    path 'versions.yml', emit: versions

    script:
    """
    longest_atg_orf.py ${transcripts} combined_transcripts.faa combined_orf_manifest.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longest_atg_orf: python-script
    END_VERSIONS
    """

    stub:
    """
    printf '>stub.1.orf1 model_id=stub.1\nMK\n' > combined_transcripts.faa
    printf 'model_id\ttranscript_id\tpeptide_id\torf_status\torf_start\torf_end\tpeptide_length\nstub.1\tstub.1\tstub.1.orf1\tCOMPLETE_ORF\t1\t9\t2\n' > combined_orf_manifest.tsv
    printf '"%s":\n    longest_atg_orf: stub\n' '${task.process}' > versions.yml
    """
}
