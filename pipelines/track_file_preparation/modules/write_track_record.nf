process WRITE_TRACK_RECORD {
    tag { "${meta.gca_accession}:${meta.id}:${meta.current_track}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${filename}"
    }

    input:
    tuple val(meta), path(track, stageAs: 'source_track'), path(checksum)

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
    cp ${track} ${artifact}
    sha=\$(awk '{print \$1}' ${checksum})
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n' > ${meta.safe_id}.${meta.current_track}.track_result.tsv
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tcomplete\t%s\ttracked\t%s\n' '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' '${meta.entity_type}' '${meta.current_track}' '${destination}' '${meta.source_path}' "\$sha" '${normalization}' >> ${meta.safe_id}.${meta.current_track}.track_result.tsv
    cp ${meta.safe_id}.${meta.current_track}.track_result.tsv ${meta.safe_id}.${meta.current_track}.provenance.tsv
    : > ${meta.safe_id}.${meta.current_track}.exceptions.tsv
    printf '"%s":\\n    record_writer: pipeline\\n' '${task.process}' > ${meta.safe_id}.${meta.current_track}.record.versions.yml
    """

    stub:
    def artifact = meta.current_track == 'coverage_forward' ? "${meta.safe_id}.forward.bw" : (meta.current_track == 'coverage_reverse' ? "${meta.safe_id}.reverse.bw" : "${meta.safe_id}.bb")
    def destination = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${artifact}"
    def normalization = meta.current_track in ['coverage_forward', 'coverage_reverse'] ?
        groovy.json.JsonOutput.toJson([tool: 'deepTools bamCoverage', bin_size: meta.coverage_bin_size, normalization: meta.coverage_normalization, strand: meta.coverage_strand]) :
        'not_applicable'
    """
    cp ${track} ${artifact}
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n${meta.gca_accession}\t${meta.assembly_release}\t${meta.id}\t${meta.entity_type}\t${meta.current_track}\t${destination}\t${meta.source_path}\tcomplete\tstub\ttracked\t${normalization}\n' > ${meta.safe_id}.${meta.current_track}.track_result.tsv
    cp ${meta.safe_id}.${meta.current_track}.track_result.tsv ${meta.safe_id}.${meta.current_track}.provenance.tsv
    : > ${meta.safe_id}.${meta.current_track}.exceptions.tsv
    printf '"%s":\\n    record_writer: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.current_track}.record.versions.yml
    """
}
