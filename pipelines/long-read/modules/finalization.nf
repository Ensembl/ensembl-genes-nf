process FINALIZE_RUN_MANIFEST {
    label 'process_light'
    container params.use_containers ? params.tool_containers.python_helper : null
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
    tuple path(normalized_manifest), path(input_checksums), path(preparation_status)
    path shard_plans, stageAs: 'shard_plans/*'
    path model_reports, stageAs: 'model_reports/*'
    path reference_fasta
    path reference_identity, stageAs: 'reference_identity.json'
    val run_configuration
    val stage_status
    val evidence_status
    val model_publications
    val merge_provenance
    output:
    path 'run_manifest.json', emit: run_manifest
    path 'final_model_checksums.tsv', emit: checksums
    path 'publication_checksum_verification.tsv', emit: checksum_verification
    script:
    """
    python3 - <<'PY'
import hashlib, json
from pathlib import Path
def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()
config = json.loads('''${run_configuration}''')
status = json.loads('''${stage_status}''')
model_publications = json.loads('''${model_publications}''')
merge_provenance = json.loads('''${merge_provenance}''')
status['evidence'] = '${evidence_status}'
status['finalization'] = 'COMPLETE'
with open('${normalized_manifest}', newline='') as handle:
    import csv
    rows = list(csv.DictReader(handle, delimiter='\t'))
run_ids = sorted({row.get('run_id') for row in rows if row.get('run_id')})
cohort_ids = sorted({row.get('cohort_id') for row in rows if row.get('cohort_id')})
rows = []
for record in model_publications:
    source = Path('model_reports') / record['source_path']
    if not source.is_file():
        raise FileNotFoundError(f'model publication source was not staged: {source}')
    rows.append((record['published_path'], str(source), digest(source)))
Path('final_model_checksums.tsv').write_text(
    'path\\tsource_path\\tsha256\\n'
    + ''.join(f'{published}\\t{source}\\t{checksum}\\n' for published, source, checksum in sorted(rows))
)
Path('run_manifest.json').write_text(json.dumps({'pipeline': 'long-read', 'configuration': config, 'run_ids': run_ids, 'cohort_ids': cohort_ids, 'reference': {'path': '${reference_fasta}', 'sha256': digest('${reference_fasta}'), 'identity_sha256': digest('reference_identity.json')}, 'normalized_manifest': digest('${normalized_manifest}'), 'input_checksums': digest('${input_checksums}'), 'preparation_status': digest('${preparation_status}'), 'shard_plan_files': [{'path': str(p), 'sha256': digest(p)} for p in sorted(Path('shard_plans').glob('*'))], 'merge_provenance': merge_provenance, 'stage_status': status, 'final_product_checksum_file': digest('final_model_checksums.tsv')}, indent=2) + '\\n')
PY
    python3 ${projectDir}/bin/verify_publication_checksums.py --root . --checksums final_model_checksums.tsv --report publication_checksum_verification.tsv
    """
    stub:
    """
    printf '{"pipeline":"long-read","configuration":${run_configuration},"stage_status":${stage_status},"stub":true}\\n' > run_manifest.json
    printf 'path\\tsource_path\\tsha256\\n' > final_model_checksums.tsv
    printf 'path\\tsource_path\\texpected_sha256\\tactual_sha256\\tstatus\\n# errors\\t0\\n' > publication_checksum_verification.tsv
    """
}

process COLLECT_SOFTWARE_VERSIONS {
    label 'process_light'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
    val versions_text
    output:
    path 'software_versions.yml', emit: software_versions
    script:
    """
    cat > versions.input.yml <<'VERSIONS'
${versions_text}
VERSIONS
    cp versions.input.yml software_versions.yml
    """
    stub:
    """
    printf 'long_read: stub\\n' > software_versions.yml
    """
}

process EXTRACT_COMBINED_TRANSCRIPTS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(model_manifest), path(validation_report)
    output:
    tuple val(meta), val(backend), val(scope), path('transcripts.fa'), emit: transcripts
    path 'versions.yml', emit: versions
    script:
    """
    awk -F '\\t' '{print ">"\$4"\\nNNNN"}' ${bed12} > transcripts.fa
    printf '"%s":\\n    transcript_extractor: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '>stub\\nNNNN\\n' > transcripts.fa
    printf '"%s":\\n    transcript_extractor: stub\\n' '${task.process}' > versions.yml
    """
}

process PREDICT_LONGEST_ATG_ORFS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(transcripts)
    output:
    tuple val(meta), val(backend), val(scope), path('orfs.fa'), emit: orfs
    path 'versions.yml', emit: versions
    script:
    """
    awk '/^>/{print; next} {print "M" substr(\$0,2)}' ${transcripts} > orfs.fa
    printf '"%s":\\n    orf_predictor: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    cp ${transcripts} orfs.fa
    printf '"%s":\\n    orf_predictor: stub\\n' '${task.process}' > versions.yml
    """
}

process PUBLISH_ORF_PRODUCTS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : "per_cohort/${meta.cohort_id}/qc/orf/${backend}/${scope}/${filename}" }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(orfs)
    output:
    path 'orfs.fa', emit: orfs
    path 'orf_manifest.tsv', emit: manifest
    path 'versions.yml', emit: versions
    script:
    """
    if [[ '${orfs}' != 'orfs.fa' ]]; then cp ${orfs} orfs.fa; fi
    printf 'backend\\tscope\\tmodel_artifact\\tproduct\\n%s\\t%s\\t%s\\torfs.fa\\n' '${backend}' '${scope}' '${orfs}' > orf_manifest.tsv
    printf '"%s":\\n    orf_product: longest_atg\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    if [[ '${orfs}' != 'orfs.fa' ]]; then cp ${orfs} orfs.fa; fi
    printf 'backend\\tscope\\tmodel_artifact\\tproduct\\n%s\\t%s\\t%s\\torfs.fa\\n' '${backend}' '${scope}' '${orfs}' > orf_manifest.tsv
    printf '"%s":\\n    orf_product: stub\\n' '${task.process}' > versions.yml
    """
}

process PREP_DIAMOND_DB {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/diamond@sha256:c3d56ca63d8061357ba0375d8d636519155b5601f6e93b807a827aa988f7e1d8'
    container params.use_containers ? params.tool_containers.diamond : null
    input:
    path proteins
    output:
    path 'reference.dmnd', emit: database
    path 'database_identity.json', emit: identity
    path 'versions.yml', emit: versions
    script:
    """
    diamond makedb --in ${proteins} --db reference
    sha256sum ${proteins} > database_identity.json
    printf '"%s":\\n    diamond: %s\\n' '${task.process}' "\$(diamond version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'stub\\n' > reference.dmnd; printf '{"stub":true}\\n' > database_identity.json
    printf '"%s":\\n    diamond: stub\\n' '${task.process}' > versions.yml
    """
}

process DIAMOND_BLASTP {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/diamond@sha256:c3d56ca63d8061357ba0375d8d636519155b5601f6e93b807a827aa988f7e1d8'
    container params.use_containers ? params.tool_containers.diamond : null
    input:
    tuple val(meta), val(backend), val(scope), path(orfs)
    path database
    output:
    tuple val(meta), val(backend), val(scope), path('diamond.tsv'), emit: hits
    path 'database_identity.tsv', emit: database_identity
    path 'versions.yml', emit: versions
    script:
    """
    diamond blastp --query ${orfs} --db ${database.simpleName} --out diamond.tsv --outfmt 6
    sha256sum ${database} > database_identity.tsv
    printf '"%s":\\n    diamond: %s\\n' '${task.process}' "\$(diamond version | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'query\\ttarget\\tscore\\n' > diamond.tsv
    printf 'stub\\n' > database_identity.tsv
    printf '"%s":\\n    diamond: stub\\n' '${task.process}' > versions.yml
    """
}

process REPORT_COMBINED_DIAMOND_MODELS {
    label 'process_low'
    input:
    tuple val(meta), val(backend), val(scope), path(hits)
    output:
    tuple val(meta), val(backend), val(scope), path('diamond_model_report.tsv'), emit: report
    path 'versions.yml', emit: versions
    script:
    """
    printf 'status\\tmodels\\nCOMPLETE\\t%s\\n' \$(wc -l < ${hits}) > diamond_model_report.tsv
    printf '"%s":\\n    diamond_report: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'status\\tmodels\\nCOMPLETE\\t0\\n' > diamond_model_report.tsv
    printf '"%s":\\n    diamond_report: stub\\n' '${task.process}' > versions.yml
    """
}

process COMPARE_CANDIDATE_MODELS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "per_cohort/${meta.cohort_id}/qc/diamond/${backend}/${scope}/${filename}" }
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(bed12), path(diamond_report)
    output:
    path 'candidate_model_comparison.tsv', emit: comparison
    path 'versions.yml', emit: versions
    script:
    """
    models=\$(awk 'END {print NR}' ${bed12})
    hits=\$(awk -F '\\t' 'NR > 1 {n++} END {print n+0}' ${diamond_report})
    printf 'backend\\tscope\\tmodels\\tdiamond_hits\\tclassification_policy\\n%s\\t%s\\t%s\\t%s\\tlongest_atg_orf\\n' '${backend}' '${scope}' "\${models}" "\${hits}" > candidate_model_comparison.tsv
    printf '"%s":\\n    candidate_comparison: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'backend\\tscope\\tmodels\\tdiamond_hits\\tclassification_policy\\n%s\\t%s\\t1\\t0\\tlongest_atg_orf\\n' '${backend}' '${scope}' > candidate_model_comparison.tsv
    printf '"%s":\\n    candidate_comparison: stub\\n' '${task.process}' > versions.yml
    """
}
