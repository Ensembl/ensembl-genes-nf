process ISOQUANT_ANNOTATION_FREE {
    label 'process_high_memory'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/isoquant/native/${filename}" }
    publishDir "${params.outdir}", mode: 'copy', enabled: params.isoquant_large_output && params.isoquant_large_output.size() > 0, saveAs: { filename -> filename.startsWith('native_products/') ? "per_run/${meta.run_id}/qc/isoquant/${filename - 'native_products/'}" : null }
    // audit container 'quay.io/biocontainers/isoquant@sha256:5b922f8211cfe6f438b814f95d83f6d141cc5ad4fc0796c3667784f102c57454'
    container params.use_containers ? params.tool_containers.isoquant : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    val data_type
    val analysis
    val high_memory
    val isoquant_args
    val large_output
    output:
    tuple val(meta), val('isoquant'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'input_manifest.json', emit: manifest
    path 'native_products/*', optional: true, emit: read_products
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: isoquant_args
    """
    isoquant.py --reference ${reference} --bam ${bam} --data_type '${data_type}' ${high_memory ? '--threads ' + task.cpus : ''} --output output ${args}
    python3 ${projectDir}/bin/finalize_gtf.py --source-dir output --pattern '*transcript_models.gtf' --output native.gtf --status status.tsv --backend isoquant --empty-message 'IsoQuant produced no transcript models' --require-source
    mkdir -p native_products
    python3 ${projectDir}/bin/collect_products.py --source output --destination native_products --patterns '*read_info*' '*read_assignments*' '*transcript_model_reads*' '*corrected_reads*' '*allinfo*' '*tagged.bam' '*deduplicated.bam'
    python3 ${projectDir}/bin/write_isoquant_manifest.py --reference ${reference} --bam ${bam} --data-type '${data_type}' --analysis ${analysis.join(' ')} --large-output ${large_output.join(' ')} --args '${args}' --output input_manifest.json
    printf '"%s":\\n    isoquant: %s\\n' '${task.process}' "\$(isoquant.py --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tisoquant\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "isoquant.1";\\n' > native.gtf
    printf 'isoquant\\tSUCCESS\\n' > status.tsv; printf '"%s":\\n    isoquant: stub\\n' '${task.process}' > versions.yml
    printf '{"stub":true}\\n' > input_manifest.json
    mkdir -p native_products
    if [[ '${large_output}' == *read_info* ]]; then printf 'stub\\n' > native_products/stub.read_info.tsv.gz; fi
    """
}

process ISOQUANT_REFERENCE_GUIDED {
    label 'process_high_memory'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/isoquant/native/${filename}" }
    publishDir "${params.outdir}", mode: 'copy', enabled: params.isoquant_large_output && params.isoquant_large_output.size() > 0, saveAs: { filename -> filename.startsWith('native_products/') ? "per_run/${meta.run_id}/qc/isoquant/${filename - 'native_products/'}" : null }
    // audit container 'quay.io/biocontainers/isoquant@sha256:5b922f8211cfe6f438b814f95d83f6d141cc5ad4fc0796c3667784f102c57454'
    container params.use_containers ? params.tool_containers.isoquant : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path genedb
    val data_type
    val analysis
    val high_memory
    val isoquant_args
    val large_output
    output:
    tuple val(meta), val('isoquant'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'input_manifest.json', emit: manifest
    path 'native_products/*', optional: true, emit: read_products
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: isoquant_args
    """
    isoquant.py --reference ${reference} --genedb ${genedb} --bam ${bam} --data_type '${data_type}' ${high_memory ? '--threads ' + task.cpus : ''} --output output ${args}
    python3 ${projectDir}/bin/finalize_gtf.py --source-dir output --pattern '*transcript_models.gtf' --output native.gtf --status status.tsv --backend isoquant --empty-message 'IsoQuant produced no transcript models' --require-source
    mkdir -p native_products
    python3 ${projectDir}/bin/collect_products.py --source output --destination native_products --patterns '*read_info*' '*read_assignments*' '*transcript_model_reads*' '*corrected_reads*' '*allinfo*' '*tagged.bam' '*deduplicated.bam'
    python3 ${projectDir}/bin/write_isoquant_manifest.py --reference ${reference} --genedb ${genedb} --bam ${bam} --data-type '${data_type}' --analysis ${analysis.join(' ')} --large-output ${large_output.join(' ')} --args '${args}' --output input_manifest.json
    printf '"%s":\\n    isoquant: %s\\n' '${task.process}' "\$(isoquant.py --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tisoquant\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "isoquant.1";\\n' > native.gtf
    printf 'isoquant\\tSUCCESS\\n' > status.tsv; printf '"%s":\\n    isoquant: stub\\n' '${task.process}' > versions.yml
    printf '{"stub":true}\\n' > input_manifest.json
    mkdir -p native_products
    if [[ '${large_output}' == *read_info* ]]; then printf 'stub\\n' > native_products/stub.read_info.tsv.gz; fi
    """
}

process NORMALISE_ISOQUANT_OUTPUTS {
    label 'process_low'
    tag "${meta.id}"
    input:
    tuple val(meta), val(backend), path(gtf), path(status)
    output:
    tuple val(meta), val(backend), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/normalize_gtf.py --gtf ${gtf} --status ${status} --output-gtf native.gtf --output-status status.tsv --backend isoquant --empty-message 'IsoQuant produced no transcript models'
    printf '"%s":\\n    isoquant_normaliser: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    if [[ '${gtf}' != 'native.gtf' ]]; then cp ${gtf} native.gtf; fi
    if [[ '${status}' != 'status.tsv' ]]; then cp ${status} status.tsv; fi
    printf '"%s":\\n    isoquant_normaliser: stub\\n' '${task.process}' > versions.yml
    """
}

process AUDIT_ISOQUANT_OUTPUTS {
    label 'process_low'
    tag "${meta.id}"
    input:
    tuple val(meta), val(backend), path(gtf), path(status)
    output:
    tuple val(meta), val(backend), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'isoquant.audit.tsv', emit: audit
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/audit_gtf.py --gtf ${gtf} --status ${status} --backend isoquant --output-gtf native.gtf --output-status status.tsv --audit isoquant.audit.tsv
    printf '"%s":\\n    isoquant_audit: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    if [[ '${gtf}' != 'native.gtf' ]]; then cp ${gtf} native.gtf; fi
    if [[ '${status}' != 'status.tsv' ]]; then cp ${status} status.tsv; fi
    printf 'backend\\tstatus\\nisoquant\\tSUCCESS\\n' > isoquant.audit.tsv
    printf '"%s":\\n    isoquant_audit: stub\\n' '${task.process}' > versions.yml
    """
}

process FLAIR_JUNCTIONS {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    input:
    tuple val(meta), val(shard_meta), path(alignment_sam)
    output:
    tuple val(meta), val(shard_meta), path('junctions.bed'), emit: junctions
    tuple val(meta), val(shard_meta), path('junction_status.tsv'), emit: status
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/bam_to_junction_bed.py --sam ${alignment_sam} --output junctions.bed
    state=HAS_JUNCTION
    if [[ ! -s junctions.bed ]]; then state=EMPTY; fi
    printf 'run_id\\tshard_id\\tstatus\\n%s\\t%s\\t%s\\n' '${meta.run_id}' '${shard_meta.shard_id}' "\${state}" > junction_status.tsv
    printf '"%s":\\n    junction_extractor: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    touch junctions.bed
    printf 'run_id\\tshard_id\\tstatus\\n%s\\t%s\\tEMPTY\\n' '${meta.run_id}' '${shard_meta.shard_id}' > junction_status.tsv
    printf '"%s":\\n    junction_extractor: stub\\n' '${task.process}' > versions.yml
    """
}

process FLAIR_TRANSCRIPTOME {
    label 'process_high_memory'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/flair/native/shards/${shard_meta.shard_id}/${filename}" }
    // audit container 'quay.io/biocontainers/flair@sha256:a39894c1507c68241fba62ddbe1a931d5a6bfaacbcb6be24631d46b60e67ce84'
    container params.use_containers ? params.tool_containers.flair : null
    input:
    tuple val(meta), val(shard_meta), path(bam), path(bai), path(junctions)
    path reference
    val flair_threads
    val flair_args
    output:
    tuple val(meta), val(shard_meta), val('flair'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'native_products/*', optional: true, emit: native_products
    path 'input_manifest.json', emit: manifest
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: flair_args
    """
    flair transcriptome -g ${reference} --sample_name '${meta.run_id}_${shard_meta.shard_id}' --junction_bed ${junctions} --junction_support 1 -b ${bam} -o flair --threads ${flair_threads ?: task.cpus} ${args}
    python3 ${projectDir}/bin/finalize_gtf.py --source flair.isoforms.gtf --output native.gtf --status status.tsv --backend flair --empty-message 'FLAIR produced no transcript models'
    mkdir -p native_products
    python3 ${projectDir}/bin/collect_products.py --source . --destination native_products --patterns 'flair*'
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products native_products --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}' --status "\${state}" --args '${args}' --input ${reference} ${bam} ${junctions}
    printf '"%s":\\n    flair: %s\\n' '${task.process}' "\$(flair --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tflair\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "flair.1";\\n' > native.gtf
    printf 'flair\\tSUCCESS\\n' > status.tsv
    printf 'chr1\\t1\\t100\\tflair.1\\t0\\t+\\n' > flair.isoforms.bed
    printf '>flair.1\\nACGT\\n' > flair.isoforms.fa
    printf 'read1\\tflair.1\\n' > flair.isoforms.read.map.txt
    mkdir -p native_products
    for product in flair*; do test -f "\${product}" && cp "\${product}" native_products/; done
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products native_products --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}' --status SUCCESS --args '${args}'
    printf '"%s":\\n    flair: stub\\n' '${task.process}' > versions.yml
    """
}

process FLAIR_EMPTY {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/flair/native/shards/${shard_meta.shard_id}/${filename}" }
    input:
    tuple val(meta), val(shard_meta), path(junction_status)
    output:
    tuple val(meta), val(shard_meta), val('flair'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'native_products/*', optional: true, emit: native_products
    path 'input_manifest.json', emit: manifest
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/require_status.py --status ${junction_status} --value EMPTY
    printf '# EMPTY: no splice junctions for shard %s\\n' '${shard_meta.shard_id}' > native.gtf
    printf 'flair\\tEMPTY\\t%s\\n' '${shard_meta.shard_id}' > status.tsv
    mkdir -p native_products
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products native_products --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}' --status EMPTY --input ${junction_status}
    printf '"%s":\\n    flair: EMPTY\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '# EMPTY: no splice junctions for shard %s\\n' '${shard_meta.shard_id}' > native.gtf
    printf 'flair\\tEMPTY\\t%s\\n' '${shard_meta.shard_id}' > status.tsv
    mkdir -p native_products
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products native_products --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}' --status EMPTY
    printf '"%s":\\n    flair: EMPTY\\n' '${task.process}' > versions.yml
    """
}

process FLAIR_COMBINE {
    label 'process_high_memory'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/flair/native/merged/${filename}" }
    input:
    tuple val(meta), val(backend), path(gtfs)
    output:
    tuple val(meta), val('flair'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'native_products/*', optional: true, emit: native_products
    path 'input_manifest.json', emit: manifest
    path 'versions.yml', emit: versions
    script:
    """
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${gtfs} --output native.gtf --backend flair --status status.tsv
    mkdir -p native_products
    python3 ${projectDir}/bin/copy_products.py --inputs ${gtfs} --destination native_products
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products native_products --run-id '${meta.run_id}' --status SUCCESS --status-file status.tsv --input ${gtfs}
    printf '"%s":\\n    flair_combine: pipeline\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${gtfs} --output native.gtf --backend flair --status status.tsv
    mkdir -p native_products
    python3 ${projectDir}/bin/copy_products.py --inputs ${gtfs} --destination native_products
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products native_products --run-id '${meta.run_id}' --status SUCCESS
    printf '"%s":\\n    flair_combine: stub\\n' '${task.process}' > versions.yml
    """
}

process BAMBU_DISCOVERY {
    label 'process_high_memory'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/bambu/native/${filename}" }
    // audit container 'quay.io/biocontainers/bioconductor-bambu@sha256:90f2d83ffb099b90bcc5a1423b2c3d2455682d96e53412a0450bfb156c4409d6'
    container params.use_containers ? params.tool_containers.bambu : null
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    val track_reads
    val quantify
    output:
    tuple val(meta), val('bambu'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'bambu.rds', emit: rds
    path 'read_assignments.tsv', emit: assignments
    path 'input_manifest.json', emit: manifest
    path 'versions.yml', emit: versions
    script:
    """
    Rscript ${projectDir}/bin/run_bambu.R ${bam} ${reference} native.gtf bambu.rds read_assignments.tsv ${track_reads} ${quantify} input_manifest.json
    python3 ${projectDir}/bin/gtf_status.py --gtf native.gtf --backend bambu --status status.tsv --empty-message 'Bambu produced no transcript models'
    printf '"%s":\\n    bambu: pinned-runtime\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'chr1\\tbambu\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "bambu.1";\\n' > native.gtf
    printf 'bambu\\tSUCCESS\\n' > status.tsv; printf '{}' > bambu.rds; printf 'read\\ttranscript\\n' > read_assignments.tsv; printf '{"backend":"bambu","stub":true}\\n' > input_manifest.json
    printf '"%s":\\n    bambu: stub\\n' '${task.process}' > versions.yml
    """
}
