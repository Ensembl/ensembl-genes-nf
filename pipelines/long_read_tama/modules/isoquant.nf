process ISOQUANT_ANNOTATION_FREE {
    tag "${meta.id}:${scope}:isoquant-annotation-free"
    label 'process_high_memory'
    container 'docker://quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0'

    input:
    tuple val(meta), path(bams), path(bais), val(scope), val(data_type), val(analysis), val(large_output), val(read_group), val(extra_args), val(report_novel_unspliced), val(check_canonical)
    path reference

    output:
    tuple val(meta), val(scope), path('isoquant_output'), path('isoquant_input_manifest.tsv'), emit: isoquant_native
    path 'versions.yml', emit: versions

    script:
    def bam_files = bams instanceof List ? bams.sort { bam_file -> bam_file.name } : [bams]
    def bam_args = bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')
    def analyses = analysis instanceof List ? analysis.join(' ') : analysis.toString()
    def outputs = large_output instanceof List ? large_output.join(' ') : large_output.toString()
    def group = read_group instanceof List ? read_group.join(' ') : read_group.toString()
    def novel = report_novel_unspliced ? '--report_novel_unspliced true' : ''
    def canonical = check_canonical ? '--check_canonical' : ''
    def numba_env = params.isoquant_numba_disable_jit ? 'NUMBA_DISABLE_JIT=1 ' : ''
    """
    mkdir -p isoquant_home
    export HOME="\$PWD/isoquant_home"
    for bam in ${bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')}; do
        test -s "\$bam" || { echo "Missing or empty IsoQuant BAM: \$bam" >&2; exit 1; }
        test -s "\${bam}.bai" || { echo "Missing or empty IsoQuant BAI for: \$bam" >&2; exit 1; }
    done
    ${numba_env}isoquant --reference '${reference}' --bam ${bam_args} --data_type '${data_type}' \\
        --analysis ${analyses} --large_output ${outputs} --read_group ${group} \\
        --output isoquant_output --prefix '${meta.id}_${scope}' --threads ${task.cpus} ${novel} ${canonical} ${extra_args}
    test -s isoquant_output/isoquant.log || { echo 'IsoQuant did not produce isoquant.log' >&2; exit 1; }
    {
        printf '%s\n' 'key\tvalue' \\
            'scope\t${scope}' \\
            'accession_or_cohort\t${meta.id}' \\
            "isoquant_version\t\$(isoquant --version)" \\
            'data_type\t${data_type}' \\
            'mode\tannotation_free' \\
            "reference_fasta_sha256\t\$(sha256sum '${reference}' | cut -d ' ' -f1)" \\
            'genedb_sha256_or_NONE\tNONE' \\
            'input_bam_names\t${bam_files.collect { bam_file -> bam_file.name }.join(',')}' \\
            "input_bam_checksums\t\$(sha256sum ${bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')} | tr '\\n' ';')"
    } > isoquant_input_manifest.tsv
    printf '"%s":\n    isoquant: \$(isoquant --version)\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p isoquant_output
    printf '# stub\nchrStub\tIsoQuant\texon\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > isoquant_output/${meta.id}_${scope}.transcript_models.gtf
    printf 'read\ttranscript\nread1\tt1\n' | gzip -c > isoquant_output/${meta.id}_${scope}.transcript_model_reads.tsv.gz
    printf 'read\tassignment\nread1\tunique\n' | gzip -c > isoquant_output/${meta.id}_${scope}.read_info.tsv.gz
    printf 'stub\n' > isoquant_output/isoquant.log
    printf 'key\tvalue\nscope\t${scope}\naccession_or_cohort\t${meta.id}\nisoquant_version\t4.0.0\ndata_type\t${data_type}\nmode\tannotation_free\nreference_fasta_sha256\tstub\ngenedb_sha256_or_NONE\tNONE\ninput_bam_names\tstub\ninput_bam_checksums\tstub\n' > isoquant_input_manifest.tsv
    printf '"%s":\n    isoquant: 4.0.0\n' '${task.process}' > versions.yml
    """
}

process ISOQUANT_REFERENCE_GUIDED {
    tag "${meta.id}:${scope}:isoquant-reference-guided"
    label 'process_high_memory'
    container 'docker://quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0'

    input:
    tuple val(meta), path(bams), path(bais), val(scope), val(data_type), val(analysis), val(large_output), val(read_group), val(extra_args), val(report_novel_unspliced), val(check_canonical)
    path reference, stageAs: 'reference.fa'
    path genedb, stageAs: 'annotation.gtf'

    output:
    tuple val(meta), val(scope), path('isoquant_output'), path('isoquant_input_manifest.tsv'), emit: isoquant_native
    path 'versions.yml', emit: versions

    script:
    def bam_files = bams instanceof List ? bams.sort { bam_file -> bam_file.name } : [bams]
    def bam_args = bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')
    def analyses = analysis instanceof List ? analysis.join(' ') : analysis.toString()
    def outputs = large_output instanceof List ? large_output.join(' ') : large_output.toString()
    def group = read_group instanceof List ? read_group.join(' ') : read_group.toString()
    def novel = report_novel_unspliced ? '--report_novel_unspliced true' : ''
    def canonical = check_canonical ? '--check_canonical' : ''
    def numba_env = params.isoquant_numba_disable_jit ? 'NUMBA_DISABLE_JIT=1 ' : ''
    def complete = params.isoquant_complete_genedb ? '--complete_genedb' : ''
    def genedb_output = params.isoquant_genedb_output ? "--genedb_output '${params.isoquant_genedb_output}'" : ''
    """
    mkdir -p isoquant_home
    export HOME="\$PWD/isoquant_home"
    for bam in ${bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')}; do
        test -s "\$bam" || { echo "Missing or empty IsoQuant BAM: \$bam" >&2; exit 1; }
        test -s "\${bam}.bai" || { echo "Missing or empty IsoQuant BAI for: \$bam" >&2; exit 1; }
    done
    ${numba_env}isoquant --reference '${reference}' --bam ${bam_args} --data_type '${data_type}' \\
        --genedb '${genedb}' ${complete} ${genedb_output} \\
        --analysis ${analyses} --large_output ${outputs} --read_group ${group} \\
        --output isoquant_output --prefix '${meta.id}_${scope}' --threads ${task.cpus} ${novel} ${canonical} ${extra_args}
    test -s isoquant_output/isoquant.log || { echo 'IsoQuant did not produce isoquant.log' >&2; exit 1; }
    {
        printf '%s\n' 'key\tvalue' \\
            'scope\t${scope}' \\
            'accession_or_cohort\t${meta.id}' \\
            "isoquant_version\t\$(isoquant --version)" \\
            'data_type\t${data_type}' \\
            'mode\treference_guided' \\
            "reference_fasta_sha256\t\$(sha256sum '${reference}' | cut -d ' ' -f1)" \\
            "genedb_sha256_or_NONE\t\$(sha256sum '${genedb}' | cut -d ' ' -f1)" \\
            'input_bam_names\t${bam_files.collect { bam_file -> bam_file.name }.join(',')}' \\
            "input_bam_checksums\t\$(sha256sum ${bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')} | tr '\\n' ';')"
    } > isoquant_input_manifest.tsv
    printf '"%s":\n    isoquant: \$(isoquant --version)\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p isoquant_output
    printf '# stub\nchrStub\tIsoQuant\texon\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > isoquant_output/${meta.id}_${scope}.transcript_models.gtf
    printf 'read\ttranscript\nread1\tt1\n' | gzip -c > isoquant_output/${meta.id}_${scope}.transcript_model_reads.tsv.gz
    printf 'read\tassignment\nread1\tunique\n' | gzip -c > isoquant_output/${meta.id}_${scope}.read_info.tsv.gz
    printf 'stub\n' > isoquant_output/isoquant.log
    printf 'key\tvalue\nscope\t${scope}\naccession_or_cohort\t${meta.id}\nisoquant_version\t4.0.0\ndata_type\t${data_type}\nmode\treference_guided\nreference_fasta_sha256\tstub\ngenedb_sha256_or_NONE\tstub\ninput_bam_names\tstub\ninput_bam_checksums\tstub\n' > isoquant_input_manifest.tsv
    printf '"%s":\n    isoquant: 4.0.0\n' '${task.process}' > versions.yml
    """
}
