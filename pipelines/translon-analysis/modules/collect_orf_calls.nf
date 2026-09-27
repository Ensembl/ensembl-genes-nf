process COLLECT_ORF_CALLS {
    label 'process_medium'
    publishDir "${params.outdir}/consensus_outputs", mode: 'copy'

    input:
    tuple val(meta), path(standardized_tsvs), path(bed12s)
    path genome_fasta
    val min_caller_agreement

    output:
    tuple val(meta), path("${meta.id}.candidate_translons.tsv"), emit: candidates
    tuple val(meta), path("${meta.id}.candidate_translons.bed12"), emit: bed12
    tuple val(meta), path("${meta.id}.characterisation_intervals.tsv"), emit: intervals
    tuple val(meta), path("${meta.id}.translation_verdicts.tsv"), emit: verdicts
    tuple val(meta), path("${meta.id}.caller_inputs/**"), emit: caller_inputs

    script:
    """
    mkdir -p ${meta.id}.caller_inputs/standardized ${meta.id}.caller_inputs/bed12
    for f in ${standardized_tsvs}; do cp "\$f" ${meta.id}.caller_inputs/standardized/; done
    for f in ${bed12s}; do cp "\$f" ${meta.id}.caller_inputs/bed12/; done
    consensus_translons.py \\
        --input-dir ${meta.id}.caller_inputs/standardized \\
        --genome-fasta ${genome_fasta} \\
        --min-caller-agreement ${min_caller_agreement} \\
        --candidates ${meta.id}.candidate_translons.tsv \\
        --bed12 ${meta.id}.candidate_translons.bed12 \\
        --intervals ${meta.id}.characterisation_intervals.tsv \\
        --verdicts ${meta.id}.translation_verdicts.tsv
    """

    stub:
    """
    mkdir -p ${meta.id}.caller_inputs
    printf 'interval_id\tchrom\tstart\tend\tstrand\tframe\tcaller_count\tcallers\tstatus\n' > ${meta.id}.candidate_translons.tsv
    printf 'stub|TX1\tchr1\t100\t220\t+\t0\t2\tribocode,ribotricer\ttrusted\n' >> ${meta.id}.candidate_translons.tsv
    printf 'chr1\t100\t220\tstub|TX1\t0\t+\t120\t0,0\t0\t120,\t0,\n' > ${meta.id}.candidate_translons.bed12
    printf 'interval_id\tchrom\tstart\tend\tstrand\tframe\nchr1:100-220:0\tchr1\t100\t220\t+\t0\n' > ${meta.id}.characterisation_intervals.tsv
    printf 'interval_id\tframe\tmechanism_class\tconfidence\tcondition_state\tpeptide_sequence\nchr1:100-220:0\t0\tcaller_consensus\t2\tcondition-unresolved\tMXX*\n' > ${meta.id}.translation_verdicts.tsv
    touch ${meta.id}.caller_inputs/stub.tsv
    """
}
