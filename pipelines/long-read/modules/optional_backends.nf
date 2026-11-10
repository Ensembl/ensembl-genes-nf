process ISOQUANT_ANNOTATION_FREE {
    label 'process_high_memory'
    // audit container 'quay.io/biocontainers/isoquant@sha256:5b922f8211cfe6f438b814f95d83f6d141cc5ad4fc0796c3667784f102c57454'
    container params.use_containers ? params.tool_containers.isoquant : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai), val(data_type)
    path reference
    val high_memory
    val isoquant_args
    output:
    tuple val(meta), path('output'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: isoquant_args
    """
    isoquant.py --reference ${reference} --bam ${bam} --data_type '${data_type}' ${high_memory ? '--threads ' + task.cpus : ''} --output output ${args}
    printf '"%s":\\n    isoquant: %s\\n' '${task.process}' "\$(isoquant.py --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    mkdir -p output
    printf 'chr1\\tisoquant\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "isoquant.1";\\n' > output/transcript_models.gtf
    touch versions.yml
    """
}

process ISOQUANT_REFERENCE_GUIDED {
    label 'process_high_memory'
    // audit container 'quay.io/biocontainers/isoquant@sha256:5b922f8211cfe6f438b814f95d83f6d141cc5ad4fc0796c3667784f102c57454'
    container params.use_containers ? params.tool_containers.isoquant : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai), val(data_type)
    path reference
    path genedb
    val high_memory
    val isoquant_args
    output:
    tuple val(meta), path('output'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: isoquant_args
    """
    isoquant.py --reference ${reference} --genedb ${genedb} --bam ${bam} --data_type '${data_type}' ${high_memory ? '--threads ' + task.cpus : ''} --output output ${args}
    printf '"%s":\\n    isoquant: %s\\n' '${task.process}' "\$(isoquant.py --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    mkdir -p output
    printf 'chr1\\tisoquant\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "isoquant.1";\\n' > output/transcript_models.gtf
    touch versions.yml
    """
}

process FINALIZE_ISOQUANT_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(isoquant_output)
    output:
    tuple val(meta), val('isoquant'), path('native.gtf'), path('status.tsv'), emit: native_models
    script:
    """
    python3 ${projectDir}/bin/finalize_gtf.py --source-dir ${isoquant_output} --pattern '*transcript_models.gtf' --output native.gtf --status status.tsv --backend isoquant --empty-message 'IsoQuant produced no transcript models' --require-source
    """
    stub:
    """
    printf 'chr1\\tisoquant\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "isoquant.1";\\n' > native.gtf
    printf 'isoquant\\tSUCCESS\\n' > status.tsv
    """
}

process COLLECT_ISOQUANT_PRODUCTS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(isoquant_output)
    val large_output
    output:
    tuple val(meta), path('native_products'), emit: products
    script:
    """
    python3 ${projectDir}/bin/collect_products.py --source ${isoquant_output} --destination native_products --patterns '*read_info*' '*read_assignments*' '*transcript_model_reads*' '*corrected_reads*' '*allinfo*' '*tagged.bam' '*deduplicated.bam'
    """
    stub:
    """
    mkdir -p native_products
    """
}

process WRITE_ISOQUANT_MANIFEST_ANNO {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai), val(data_type)
    path reference
    val analysis
    val large_output
    val isoquant_args
    output:
    path 'input_manifest.json', emit: manifest
    script:
    """
    python3 ${projectDir}/bin/write_isoquant_manifest.py --reference ${reference} --bam ${bam} --data-type '${data_type}' --analysis ${analysis.join(' ')} --large-output ${large_output.join(' ')} --args '${isoquant_args}' --output input_manifest.json
    """
    stub:
    """
    printf '{"stub":true}\\n' > input_manifest.json
    """
}

process WRITE_ISOQUANT_MANIFEST_GUIDED {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai), val(data_type)
    path reference
    path genedb
    val analysis
    val large_output
    val isoquant_args
    output:
    path 'input_manifest.json', emit: manifest
    script:
    """
    python3 ${projectDir}/bin/write_isoquant_manifest.py --reference ${reference} --genedb ${genedb} --bam ${bam} --data-type '${data_type}' --analysis ${analysis.join(' ')} --large-output ${large_output.join(' ')} --args '${isoquant_args}' --output input_manifest.json
    """
    stub:
    """
    printf '{"stub":true}\\n' > input_manifest.json
    """
}

process NORMALISE_ISOQUANT_OUTPUTS {
    label 'process_low'
    tag "${meta.id}"
    input:
    tuple val(meta), val(backend), path(gtf), path(status)
    output:
    tuple val(meta), val(backend), path('native.gtf'), path('status.tsv'), emit: native_models
    script:
    """
    python3 ${projectDir}/bin/normalize_gtf.py --gtf ${gtf} --status ${status} --output-gtf native.gtf --output-status status.tsv --backend isoquant --empty-message 'IsoQuant produced no transcript models'
    """
    stub:
    """
    if [[ "${gtf}" != 'native.gtf' ]]; then cp ${gtf} native.gtf; fi
    if [[ "${status}" != 'status.tsv' ]]; then cp ${status} status.tsv; fi
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
    script:
    """
    python3 ${projectDir}/bin/audit_gtf.py --gtf ${gtf} --status ${status} --backend isoquant --output-gtf native.gtf --output-status status.tsv --audit isoquant.audit.tsv
    """
    stub:
    """
    if [[ "${gtf}" != 'native.gtf' ]]; then cp ${gtf} native.gtf; fi
    if [[ "${status}" != 'status.tsv' ]]; then cp ${status} status.tsv; fi
    printf 'backend\\tstatus\\nisoquant\\tSUCCESS\\n' > isoquant.audit.tsv
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
    script:
    """
    python3 ${projectDir}/bin/bam_to_junction_bed.py --sam ${alignment_sam} --output junctions.bed --status junction_status.tsv --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}'
    """
    stub:
    """
    touch junctions.bed
    printf 'run_id\\tshard_id\\tstatus\\n%s\\t%s\\tEMPTY\\n' '${meta.run_id}' '${shard_meta.shard_id}' > junction_status.tsv
    """
}

process FLAIR_TRANSCRIPTOME {
    label 'process_high_memory'
    // audit container 'quay.io/biocontainers/flair@sha256:a39894c1507c68241fba62ddbe1a931d5a6bfaacbcb6be24631d46b60e67ce84'
    container params.use_containers ? params.tool_containers.flair : null
    input:
    tuple val(meta), val(shard_meta), path(bam), path(bai), path(junctions)
    path reference
    val flair_threads
    val flair_args
    output:
    tuple val(meta), val(shard_meta), path('flair*'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: flair_args
    """
    flair transcriptome -g ${reference} --sample_name '${meta.run_id}_${shard_meta.shard_id}' --junction_bed ${junctions} --junction_support 1 -b ${bam} -o flair --threads ${flair_threads ?: task.cpus} ${args}
    printf '"%s":\\n    flair: %s\\n' '${task.process}' "\$(flair --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tflair\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "flair.1";\\n' > flair.isoforms.gtf
    printf 'chr1\\t1\\t100\\tflair.1\\t0\\t+\\n' > flair.isoforms.bed
    printf '>flair.1\\nACGT\\n' > flair.isoforms.fa
    printf 'read1\\tflair.1\\n' > flair.isoforms.read.map.txt
    touch versions.yml
    """
}

process FINALIZE_FLAIR_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(flair_products)
    output:
    tuple val(meta), val(shard_meta), val('flair'), path('*.flair.gtf'), path('*.flair.status.tsv'), emit: native_models
    script:
    def flair_gtf = "${meta.run_id}.${shard_meta.shard_id}.flair.gtf"
    def flair_status = "${meta.run_id}.${shard_meta.shard_id}.flair.status.tsv"
    """
    python3 ${projectDir}/bin/finalize_gtf.py --source flair.isoforms.gtf --output '${flair_gtf}' --status '${flair_status}' --backend flair --empty-message 'FLAIR produced no transcript models'
    """
    stub:
    def flair_gtf = "${meta.run_id}.${shard_meta.shard_id}.flair.gtf"
    def flair_status = "${meta.run_id}.${shard_meta.shard_id}.flair.status.tsv"
    """
    printf 'chr1\\tflair\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "flair.1";\\n' > '${flair_gtf}'
    printf 'flair\\tSUCCESS\\n' > '${flair_status}'
    """
}

process COLLECT_FLAIR_PRODUCTS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(flair_products)
    output:
    tuple val(meta), val(shard_meta), path('native_products'), emit: products
    script:
    """
    python3 ${projectDir}/bin/collect_products.py --source . --destination native_products --patterns 'flair*'
    """
    stub:
    """
    mkdir -p native_products
    printf 'stub\\n' > native_products/flair.stub
    """
}

process WRITE_FLAIR_MANIFEST {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(native_products), path(status), path(bam), path(junctions)
    path reference
    val flair_args
    output:
    path 'input_manifest.json', emit: manifest
    script:
    """
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products ${native_products} --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}' --status SUCCESS --status-file ${status} --args '${flair_args}' --input ${reference} ${bam} ${junctions}
    """
    stub:
    """
    printf '{"backend":"flair","stub":true}\\n' > input_manifest.json
    """
}

process FLAIR_EMPTY {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    input:
    tuple val(meta), val(shard_meta), path(junction_status)
    output:
    tuple val(meta), val(shard_meta), val('flair'), path('*.flair.gtf'), path('*.flair.status.tsv'), emit: native_models
    path 'native_products', emit: native_products
    path 'input_manifest.json', emit: manifest
    script:
    def flair_gtf = "${meta.run_id}.${shard_meta.shard_id}.flair.gtf"
    def flair_status = "${meta.run_id}.${shard_meta.shard_id}.flair.status.tsv"
    """
    python3 ${projectDir}/bin/write_empty_flair.py --junction-status ${junction_status} --run-id '${meta.run_id}' --shard-id '${shard_meta.shard_id}' --gtf '${flair_gtf}' --status '${flair_status}' --manifest input_manifest.json --products native_products
    """
    stub:
    def flair_gtf = "${meta.run_id}.${shard_meta.shard_id}.flair.gtf"
    def flair_status = "${meta.run_id}.${shard_meta.shard_id}.flair.status.tsv"
    """
    printf '# EMPTY: no splice junctions for shard %s\\n' '${shard_meta.shard_id}' > '${flair_gtf}'
    printf 'flair\\tEMPTY\\t%s\\n' '${shard_meta.shard_id}' > '${flair_status}'
    mkdir -p native_products
    printf '{"backend":"flair","status":"EMPTY"}\\n' > input_manifest.json
    """
}

process FLAIR_COMBINE_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), val(backend), path(gtfs)
    output:
    tuple val(meta), val('flair'), path('native.gtf'), path('status.tsv'), emit: native_models
    script:
    """
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${gtfs} --output native.gtf --backend flair --status status.tsv
    """
    stub:
    """
    printf 'chr1\\tflair\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "flair.1";\\n' > native.gtf
    printf 'backend\\tstatus\\nflair\\tSUCCESS\\n' > status.tsv
    """
}

process FLAIR_COMBINE_PRODUCTS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), val(backend), path(gtfs)
    output:
    tuple val(meta), path(gtfs), path('native_products'), emit: products
    script:
    """
    python3 ${projectDir}/bin/copy_products.py --inputs ${gtfs} --destination native_products
    """
    stub:
    """
    mkdir -p native_products
    printf 'stub\\n' > native_products/flair.stub
    """
}

process WRITE_FLAIR_COMBINED_MANIFEST {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(native_products), path(status), path(gtfs)
    output:
    path 'input_manifest.json', emit: manifest
    script:
    """
    python3 ${projectDir}/bin/write_flair_manifest.py --output input_manifest.json --native-products ${native_products} --run-id '${meta.run_id}' --status SUCCESS --status-file ${status} --input ${gtfs}
    """
    stub:
    """
    printf '{"backend":"flair","stub":true}\\n' > input_manifest.json
    """
}

process BAMBU_DISCOVERY {
    label 'process_high_memory'
    // audit container 'quay.io/biocontainers/bioconductor-bambu@sha256:90f2d83ffb099b90bcc5a1423b2c3d2455682d96e53412a0450bfb156c4409d6'
    container params.use_containers ? params.tool_containers.bambu : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    val track_reads
    val quantify
    output:
    tuple val(meta), path('native.gtf'), emit: raw
    path 'bambu.rds', emit: rds
    path 'read_assignments.tsv', emit: assignments
    path 'input_manifest.json', emit: manifest
    script:
    """
    Rscript ${projectDir}/bin/run_bambu.R ${bam} ${reference} native.gtf bambu.rds read_assignments.tsv ${track_reads} ${quantify} input_manifest.json
    """
    stub:
    """
    printf 'chr1\\tbambu\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "bambu.1";\\n' > native.gtf
    printf '{}' > bambu.rds
    printf 'read\\ttranscript\\n' > read_assignments.tsv
    printf '{"backend":"bambu","stub":true}\\n' > input_manifest.json
    """
}

process BAMBU_STATUS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(native_gtf)
    output:
    tuple val(meta), val('bambu'), path('native.gtf'), path('status.tsv'), emit: native_models
    script:
    """
    python3 ${projectDir}/bin/gtf_status.py --gtf ${native_gtf} --backend bambu --status status.tsv --empty-message 'Bambu produced no transcript models'
    """
    stub:
    """
    if [[ "${native_gtf}" != 'native.gtf' ]]; then cp ${native_gtf} native.gtf; fi
    printf 'bambu\\tSUCCESS\\n' > status.tsv
    """
}

process TRANSGRAM_ANNOTATION_FREE {
    label 'process_high_memory'
    // Official project-provided image: trans_gram-20250912.sif
    container params.use_containers ? params.tool_containers.transgram : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai), val(data_type)
    val transgram_args
    val transgram_executable
    output:
    tuple val(meta), path('transgram_outdir'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: transgram_args
    def data_type_arg = data_type == 'auto' ? '' : "--specify_data_type ${data_type}"
    def executable = task.ext.executable ?: transgram_executable
    """
    ${executable} --bam ${bam} --output_dir transgram_outdir ${data_type_arg} ${args}
    test -f transgram_outdir/TransGram.gtf
    test -f transgram_outdir/TransGram-singleExon.gtf
    version_status=0
    version_output=\$(${executable} --version 2>&1) || version_status=\$?
    if [[ \${version_status} -ne 0 && \${version_status} -ne 1 ]]; then
        echo "TransGram version check failed with exit status \${version_status}" >&2
        exit \${version_status}
    fi
    printf '"%s":\\n    transgram: %s\\n' '${task.process}' "\${version_output}" > versions.yml
    """
    stub:
    """
    mkdir -p transgram_outdir
    printf 'chr1\\tTransGram\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "transgram.1";\\n' > transgram_outdir/TransGram.gtf
    printf 'chr1\\tTransGram\\texon\\t200\\t250\\t.\\t+\\t.\\tgene_id "g2"; transcript_id "transgram.single.1";\\n' > transgram_outdir/TransGram-singleExon.gtf
    printf 'transcript\\tabundance\\ntransgram.1\\t1\\n' > transgram_outdir/TransGram-Abundance.tsv
    touch versions.yml
    """
}

process TRANSGRAM_ANNOTATION_GUIDED {
    label 'process_high_memory'
    // Official project-provided image: trans_gram-20250912.sif
    container params.use_containers ? params.tool_containers.transgram : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(bam), path(bai), val(data_type)
    path annotation
    val transgram_args
    val transgram_executable
    output:
    tuple val(meta), path('transgram_outdir'), emit: raw
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: transgram_args
    def data_type_arg = data_type == 'auto' ? '' : "--specify_data_type ${data_type}"
    def executable = task.ext.executable ?: transgram_executable
    """
    ${executable} --bam ${bam} --annotation ${annotation} --output_dir transgram_outdir ${data_type_arg} ${args}
    test -f transgram_outdir/TransGram.gtf
    test -f transgram_outdir/TransGram-singleExon.gtf
    version_status=0
    version_output=\$(${executable} --version 2>&1) || version_status=\$?
    if [[ \${version_status} -ne 0 && \${version_status} -ne 1 ]]; then
        echo "TransGram version check failed with exit status \${version_status}" >&2
        exit \${version_status}
    fi
    printf '"%s":\\n    transgram: %s\\n' '${task.process}' "\${version_output}" > versions.yml
    """
    stub:
    """
    mkdir -p transgram_outdir
    printf 'chr1\\tTransGram\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "transgram.1";\\n' > transgram_outdir/TransGram.gtf
    printf 'chr1\\tTransGram\\texon\\t200\\t250\\t.\\t+\\t.\\tgene_id "g2"; transcript_id "transgram.single.1";\\n' > transgram_outdir/TransGram-singleExon.gtf
    printf 'transcript\\tabundance\\ntransgram.1\\t1\\n' > transgram_outdir/TransGram-Abundance.tsv
    touch versions.yml
    """
}

process FINALIZE_TRANSGRAM_GTF {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}"
    input:
    tuple val(meta), path(transgram_output)
    output:
    tuple val(meta), val('transgram'), path('native.gtf'), path('status.tsv'), emit: native_models
    script:
    """
    test -f ${transgram_output}/TransGram.gtf
    test -f ${transgram_output}/TransGram-singleExon.gtf
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${transgram_output}/TransGram.gtf ${transgram_output}/TransGram-singleExon.gtf --output native.gtf --backend transgram --status status.tsv
    """
    stub:
    """
    printf 'chr1\\tTransGram\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "transgram.1";\\n' > native.gtf
    printf 'backend\\tstatus\\ntransgram\\tCOMPLETE\\n' > status.tsv
    """
}
