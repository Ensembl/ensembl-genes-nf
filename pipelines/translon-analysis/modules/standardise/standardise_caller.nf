process STANDARDISE_CALLER {
    tag "${meta.id}:${tool}:${meta.codon ?: 'all'}:${meta.shard_id ?: 'all'}"
    label 'process_low'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/ribotricer:1.5.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(rawdir)
    val tool
    path gtf

    output:
    tuple val(meta), val(tool), path("${meta.id}_${tool}.bed12"), emit: bed12
    tuple val(meta), val(tool), path("${meta.id}_${tool}.tsv"), emit: standardized
    tuple val(meta), val(tool), path("${meta.id}_${tool}.status.json"), emit: status

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${params.translon_analysis_bin}/standardise_caller.py \\
        --adapter-version 2 --raw ${rawdir} --tool ${tool} --sample ${meta.id} \\
        --gtf ${gtf} --output ${meta.id}_${tool}.tsv --bed12 ${meta.id}_${tool}.bed12 \\
        --status ${meta.id}_${tool}.status.json \\
        --requested-start-codons '${params.start_codons}' \\
        --requested-stop-codons '${params.stop_codons}' --caller-codon '${meta.codon ?: ''}' --shard-id '${meta.shard_id ?: 'all'}' ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardiser: python3
    END_VERSIONS
    """

    stub:
    """
    printf 'sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\ttranscript_id\torf_id\tscore\tpval\tqval\textra_json\n${meta.id}\t${tool}\tchr1\t100\t200\t+\t0\tTX1\tORF1\t10\t0.01\t\t{}\n' > ${meta.id}_${tool}.tsv
    printf 'chr1\t100\t200\tORF1|TX1|${tool}\t900\t+\t100\t200\t0,0,0\t1\t100,\t0,\n' > ${meta.id}_${tool}.bed12
    printf '{"sample_id":"${meta.id}","tool":"${tool}","status":"ok","standardised_rows":1}\n' > ${meta.id}_${tool}.status.json
    printf '"stub":\n    standardiser: stub\n' > versions.yml
    """
}
