process COLLECT_ORF_CALLS {
    label 'process_medium'
    publishDir "${params.outdir}/01_consensus", mode: 'copy'

    input:
    path standardized_tsvs
    path bed12s
    path genome_fasta
    val min_caller_agreement

    output:
    path 'candidate_translons.tsv', emit: candidates
    path 'candidate_translons.bed12', emit: bed12
    path 'characterisation_intervals.tsv', emit: intervals
    path 'translation_verdicts.tsv', emit: verdicts
    path 'caller_inputs/**', emit: caller_inputs

    script:
    """
    mkdir -p caller_inputs
    for f in ${standardized_tsvs}; do cp "\$f" caller_inputs/; done
    for f in ${bed12s}; do cp "\$f" caller_inputs/; done
    consensus_translons.py \\
        --input-dir caller_inputs \\
        --genome-fasta ${genome_fasta} \\
        --min-caller-agreement ${min_caller_agreement} \\
        --candidates candidate_translons.tsv \\
        --bed12 candidate_translons.bed12 \\
        --intervals characterisation_intervals.tsv \\
        --verdicts translation_verdicts.tsv
    """

    stub:
    """
    mkdir -p caller_inputs
    printf 'interval_id\tchrom\tstart\tend\tstrand\tframe\tcaller_count\tcallers\tstatus\n' > candidate_translons.tsv
    printf 'stub|TX1\tchr1\t100\t220\t+\t0\t2\tribocode,ribotricer\ttrusted\n' >> candidate_translons.tsv
    printf 'chr1\t100\t220\tstub|TX1\t0\t+\t120\t0,0\t0\t120,\t0,\n' > candidate_translons.bed12
    printf 'interval_id\tchrom\tstart\tend\tstrand\tframe\nchr1:100-220:0\tchr1\t100\t220\t+\t0\n' > characterisation_intervals.tsv
    printf 'interval_id\tframe\tmechanism_class\tconfidence\tcondition_state\tpeptide_sequence\nchr1:100-220:0\t0\tcaller_consensus\t2\tcondition-unresolved\tMXX*\n' > translation_verdicts.tsv
    touch caller_inputs/stub.tsv
    """
}
