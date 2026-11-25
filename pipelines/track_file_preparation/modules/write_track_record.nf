process WRITE_TRACK_RECORD {
    tag { "${meta.gca_accession}:${meta.id}:${meta.current_track}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${filename}"
    }

    input:
    tuple val(meta), path(track, stageAs: 'source_track'), path(checksum)
    path record_script

    output:
    tuple val(meta), path("${meta.current_track == 'coverage_forward' ? meta.safe_id + '.forward.bw' : (meta.current_track == 'coverage_reverse' ? meta.safe_id + '.reverse.bw' : meta.safe_id + '.bb')}"), path("${meta.safe_id}.${meta.current_track}.track_result.tsv"), path("${meta.safe_id}.${meta.current_track}.provenance.tsv"), path("${meta.safe_id}.${meta.current_track}.exceptions.tsv"), emit: results
    path "${meta.safe_id}.${meta.current_track}.record.versions.yml", emit: versions

    script:
    def artifact = meta.current_track == 'coverage_forward' ? "${meta.safe_id}.forward.bw" : (meta.current_track == 'coverage_reverse' ? "${meta.safe_id}.reverse.bw" : "${meta.safe_id}.bb")
    def destination = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${artifact}"
    def normalization = meta.current_track in ['coverage_forward', 'coverage_reverse'] ?
        groovy.json.JsonOutput.toJson([tool: 'deepTools bamCoverage', bin_size: meta.coverage_bin_size, normalization: meta.coverage_normalization, strand: meta.coverage_strand]) :
        'not_applicable'
    """
    python3 ${record_script} \\
        --track ${track} \\
        --checksum ${checksum} \\
        --artifact '${artifact}' \\
        --track-result '${meta.safe_id}.${meta.current_track}.track_result.tsv' \\
        --provenance '${meta.safe_id}.${meta.current_track}.provenance.tsv' \\
        --exceptions '${meta.safe_id}.${meta.current_track}.exceptions.tsv' \\
        --versions '${meta.safe_id}.${meta.current_track}.record.versions.yml' \\
        --gca '${meta.gca_accession}' \\
        --assembly-release '${meta.assembly_release}' \\
        --entity-id '${meta.id}' \\
        --entity-type '${meta.entity_type}' \\
        --track-type '${meta.current_track}' \\
        --destination '${destination}' \\
        --source-path '${meta.source_path}' \\
        --normalization '${normalization}'
    """

    stub:
    def artifact = meta.current_track == 'coverage_forward' ? "${meta.safe_id}.forward.bw" : (meta.current_track == 'coverage_reverse' ? "${meta.safe_id}.reverse.bw" : "${meta.safe_id}.bb")
    def destination = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${artifact}"
    def normalization = meta.current_track in ['coverage_forward', 'coverage_reverse'] ?
        groovy.json.JsonOutput.toJson([tool: 'deepTools bamCoverage', bin_size: meta.coverage_bin_size, normalization: meta.coverage_normalization, strand: meta.coverage_strand]) :
        'not_applicable'
    """
    python3 ${record_script} \\
        --track ${track} \\
        --checksum ${checksum} \\
        --artifact '${artifact}' \\
        --track-result '${meta.safe_id}.${meta.current_track}.track_result.tsv' \\
        --provenance '${meta.safe_id}.${meta.current_track}.provenance.tsv' \\
        --exceptions '${meta.safe_id}.${meta.current_track}.exceptions.tsv' \\
        --versions '${meta.safe_id}.${meta.current_track}.record.versions.yml' \\
        --gca '${meta.gca_accession}' \\
        --assembly-release '${meta.assembly_release}' \\
        --entity-id '${meta.id}' \\
        --entity-type '${meta.entity_type}' \\
        --track-type '${meta.current_track}' \\
        --destination '${destination}' \\
        --source-path '${meta.source_path}' \\
        --normalization '${normalization}' \\
        --stub
    """
}
