process TAMA_PREPARE_SAM {
    label 'process_low'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    output:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path('shard.sam'), emit: prepared
    path 'versions.yml', emit: versions
    script:
    """
    samtools view -h -F 3844 ${shard_bam} | awk '!/^@PG/' > shard.sam
    printf '"%s":\\n    samtools: %s\\n' '${task.process}' "\$(samtools --version | head -n1)" > versions.yml
    """
    stub:
    """
    touch shard.sam
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}

process TAMA_COLLAPSE {
    label 'process_high_memory'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/tama/native/shards/${shard_meta.shard_id}/${filename}" }
    // audit container 'docker.io/sguizard/tama@sha256:66a05355155f5ecef2f585736e7fe3d9d9ad38b566f1179437828ad2d7bfd123'
    container params.use_containers ? params.tool_containers.tama : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai), path(shard_sam)
    path reference
    val cap_mode
    val end_mode
    val run_mode
    output:
    tuple val(meta), val(shard_meta), val('tama'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'input_checksums.tsv', emit: checksums
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: "-x ${cap_mode} -e ${end_mode} -rm ${run_mode}"
    """
    tama_collapse.py -s ${shard_sam} -f ${reference} -p tama -m 10 ${args}
    if test -s tama.bed; then
        awk -f ${projectDir}/bin/bed12_to_gtf.awk tama.bed > native.gtf
    else
        printf '# EMPTY: TAMA produced no models for shard %s\\n' '${shard_meta.shard_id}' > native.gtf
    fi
    test -s native.gtf
    if grep -q '^# EMPTY:' native.gtf; then state=EMPTY; else state=SUCCESS; fi
    printf 'backend\\tstatus\\tshard_id\\n' > status.tsv
    printf 'tama\\t%s\\t%s\\n' "\${state}" '${shard_meta.shard_id}' >> status.tsv
    sha256sum ${shard_bam} ${shard_bai} ${shard_sam} > input_checksums.tsv
    cat > versions.yml <<VERSIONS
"${task.process}":
    tama: \$(tama_collapse.py -v 1 2>&1 | head -n1)
VERSIONS
    """
    stub:
    """
    printf 'chr1\\ttama\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "tama.1";\\n' > native.gtf
    printf 'backend\\tstatus\\tshard_id\\ntama\\tSUCCESS\\t%s\\n' '${shard_meta.shard_id}' > status.tsv
    printf 'stub\\n' > input_checksums.tsv
    printf '"%s":\\n    tama: stub\\n' '${task.process}' > versions.yml
    """
}

process STRINGTIE2_COLLAPSE {
    label 'process_high'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/stringtie2/native/shards/${shard_meta.shard_id}/${filename}" }
    // audit container 'quay.io/biocontainers/stringtie@sha256:84207eb013dbcff90928d8ddc75902580d3885ce8edaa36b8fab7665a3dc2e3c'
    container params.use_containers ? params.tool_containers.stringtie2 : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    path reference
    output:
    tuple val(meta), val(shard_meta), val('stringtie2'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'input_checksums.tsv', emit: checksums
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: ''
    """
    stringtie ${shard_bam} -o native.gtf ${args}
    test -s native.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' native.gtf; then state=SUCCESS; else printf '# EMPTY: StringTie2 produced no models for shard %s\\n' '${shard_meta.shard_id}' > native.gtf; state=EMPTY; fi
    printf 'stringtie2\\t%s\\t%s\\n' "\${state}" '${shard_meta.shard_id}' > status.tsv
    sha256sum ${shard_bam} ${shard_bai} > input_checksums.tsv
    printf '"%s":\\n    stringtie: %s\\n' '${task.process}' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tstringtie2\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "stringtie2.1";\\n' > native.gtf
    printf 'stringtie2\\tSUCCESS\\t%s\\n' '${shard_meta.shard_id}' > status.tsv; printf 'stub\\n' > input_checksums.tsv
    printf '"%s":\\n    stringtie: stub\\n' '${task.process}' > versions.yml
    """
}

process STRINGTIE3_COLLAPSE {
    label 'process_high'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/stringtie3/native/shards/${shard_meta.shard_id}/${filename}" }
    // audit container 'quay.io/biocontainers/stringtie@sha256:1e2d08591607f70ee10b830a5fd17eaa3c1c34d25ca28cb224a73658053a1f31'
    container params.use_containers ? params.tool_containers.stringtie3 : null
    tag "${shard_meta.shard_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    path reference
    output:
    tuple val(meta), val(shard_meta), val('stringtie3'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'input_checksums.tsv', emit: checksums
    path 'versions.yml', emit: versions
    script:
    """
    stringtie --version > /dev/null
    stringtie ${shard_bam} -o native.gtf
    test -s native.gtf
    if awk -F '\\t' '\$3 == "exon" { found=1 } END { exit(found ? 0 : 1) }' native.gtf; then state=SUCCESS; else printf '# EMPTY: StringTie3 produced no models for shard %s\\n' '${shard_meta.shard_id}' > native.gtf; state=EMPTY; fi
    printf 'stringtie3\\t%s\\t%s\\n' "\${state}" '${shard_meta.shard_id}' > status.tsv; sha256sum ${shard_bam} ${shard_bai} > input_checksums.tsv
    printf '"%s":\\n    stringtie: %s\\n' '${task.process}' "\$(stringtie --version 2>&1 | head -n1)" > versions.yml
    """
    stub:
    """
    printf 'chr1\\tstringtie3\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "stringtie3.1";\\n' > native.gtf
    printf 'stringtie3\\tSUCCESS\\t%s\\n' '${shard_meta.shard_id}' > status.tsv; printf 'stub\\n' > input_checksums.tsv
    printf '"%s":\\n    stringtie: stub\\n' '${task.process}' > versions.yml
    """
}

process BAM_TO_ALIGNMENT_GTF {
    label 'process_medium'
    // audit container 'quay.io/biocontainers/samtools@sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5'
    container params.use_containers ? params.tool_containers.samtools : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(shard_meta), path(shard_bam), path(shard_bai)
    output:
    tuple val(meta), val(shard_meta), path('alignment.gtf'), emit: alignment_gtf
    path 'versions.yml', emit: versions
    script:
    """
    samtools view -F 3844 ${shard_bam} | awk -F '\\t' '
        BEGIN { OFS = "\\t" }
        {
            read_id = \$1; flag = \$2 + 0; chrom = \$3; reference_position = \$4 - 1
            strand = (int(flag / 16) % 2 ? "-" : "+")
            exon_start = reference_position; exon_number = 0; cigar = \$6; rest = cigar
            while (match(rest, /[0-9]+[MIDNSHP=X]/)) {
                token = substr(rest, RSTART, RLENGTH)
                length_bp = token + 0; operation = substr(token, length(token), 1)
                if (operation ~ /[M=X]/) {
                    reference_position += length_bp
                } else if (operation == "N") {
                    if (reference_position > exon_start) {
                        exon_number++
                        printf "%s\\tpipeline\\texon\\t%d\\t%d\\t.\\t%s\\t.\\tgene_id \\\"%s\\\"; transcript_id \\\"%s\\\"; exon_number \\\"%d\\\";\\n", chrom, exon_start + 1, reference_position, strand, read_id, read_id, exon_number
                        written++
                    }
                    reference_position += length_bp; exon_start = reference_position
                } else if (operation == "D") {
                    reference_position += length_bp
                }
                rest = substr(rest, RSTART + RLENGTH)
            }
            if (reference_position > exon_start) {
                exon_number++
                printf "%s\\tpipeline\\texon\\t%d\\t%d\\t.\\t%s\\t.\\tgene_id \\\"%s\\\"; transcript_id \\\"%s\\\"; exon_number \\\"%d\\\";\\n", chrom, exon_start + 1, reference_position, strand, read_id, read_id, exon_number
                written++
            }
        }
        END { if (written == 0) exit 1 }
    ' > alignment.gtf
    test -s alignment.gtf
    cat > versions.yml <<VERSIONS
"${task.process}":
    converter: pipeline
VERSIONS
    """
    stub:
    """
    printf 'chr1\\tconverter\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "read.1";\\n' > alignment.gtf
    printf '"%s":\\n    converter: stub\\n' '${task.process}' > versions.yml
    """
}

process TMERGE_COLLAPSE {
    label 'process_high'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_native_intermediates, saveAs: { filename -> filename == 'versions.yml' ? null : "per_run/${meta.run_id}/models/tmerge/native/shards/${shard_meta.shard_id}/${filename}" }
    // audit container 'docker.io/library/perl@sha256:587c0c65838dca2fa760e0a0ab169f7260894f208e68f38535ca644f397470c7'
    container params.use_containers ? params.tool_containers.tmerge : null
    tag "${meta.run_id}"
    input:
    tuple val(meta), val(shard_meta), path(alignment_gtf)
    val tmerge_args
    output:
    tuple val(meta), val(shard_meta), val('tmerge'), path('native.gtf'), path('status.tsv'), emit: native_models
    path 'input_checksums.tsv', emit: checksums
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: tmerge_args
    """
    perl ${projectDir}/bin/tmerge ${args} ${alignment_gtf} > native.gtf
    if test -s native.gtf; then state=SUCCESS; else printf '# EMPTY: tmerge produced no models for shard %s\\n' '${shard_meta.shard_id}' > native.gtf; state=EMPTY; fi
    printf 'tmerge\\t%s\\t%s\\n' "\${state}" '${shard_meta.shard_id}' > status.tsv; sha256sum ${alignment_gtf} > input_checksums.tsv
    printf '"%s":\\n    tmerge: vendored-v1.0-8b4d6e7c1c94955931946081476e326c4cece161\\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf 'chr1\\ttmerge\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "tmerge.1";\\n' > native.gtf
    printf 'tmerge\\tSUCCESS\\t%s\\n' '${shard_meta.shard_id}' > status.tsv; printf 'stub\\n' > input_checksums.tsv
    printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
