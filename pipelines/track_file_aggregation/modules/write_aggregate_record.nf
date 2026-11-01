process WRITE_AGGREGATE_RECORD {
    tag { "${meta.gca_accession}:${meta.aggregate_id}:${meta.current_track}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${filename}"
    }

    input:
    tuple val(meta), path(track), path(checksum)
    path record_script

    output:
    tuple val(meta), path("${meta.output_name}"), path("${meta.safe_id}.${meta.current_track}.track_result.tsv"), path("${meta.safe_id}.${meta.current_track}.provenance.json"), emit: results
    path "${meta.safe_id}.${meta.current_track}.record.versions.yml", emit: versions

    script:
    def published_track_path = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${meta.output_name}"
    def provenance_path = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${meta.safe_id}.${meta.current_track}.provenance.json"
    def metadata_json = groovy.json.JsonOutput.toJson(meta + [published_track_path: published_track_path, provenance_path: provenance_path])
    """
    cp ${track} ${meta.output_name}
    cat > aggregate_record_metadata.json <<'JSON'
    ${metadata_json}
    JSON
    python3 ${record_script} \
        --metadata aggregate_record_metadata.json \
        --track ${meta.output_name} \
        --checksum ${checksum} \
        --result ${meta.safe_id}.${meta.current_track}.track_result.tsv \
        --provenance ${meta.safe_id}.${meta.current_track}.provenance.json
    printf '"%s":\\n    aggregate_record: pipeline\\n' '${task.process}' > ${meta.safe_id}.${meta.current_track}.record.versions.yml
    """

    stub:
    def published_track_path = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${meta.output_name}"
    def provenance_path = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/${meta.current_track}/${meta.safe_id}.${meta.current_track}.provenance.json"
    def metadata_json = groovy.json.JsonOutput.toJson(meta + [published_track_path: published_track_path, provenance_path: provenance_path])
    """
    printf 'stub\\n' > ${meta.output_name}
    printf 'stub\\n' > ${meta.safe_id}.${meta.current_track}.sha256
    cat > aggregate_record_metadata.json <<'JSON'
    ${metadata_json}
    JSON
    python3 ${record_script} \
        --metadata aggregate_record_metadata.json \
        --track ${meta.output_name} \
        --checksum ${meta.safe_id}.${meta.current_track}.sha256 \
        --result ${meta.safe_id}.${meta.current_track}.track_result.tsv \
        --provenance ${meta.safe_id}.${meta.current_track}.provenance.json
    printf '"%s":\\n    aggregate_record: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.current_track}.record.versions.yml
    """
}
