process RUN_ORFQUANT {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfquant:1.1.0-txdbmaker2'
    input:
    tuple val(meta), path(bam), path(bai), path(offsets)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    export MPLCONFIGDIR=\$PWD/.mplconfig
    mkdir -p raw \"\$MPLCONFIGDIR\"
    Rscript ${params.translon_analysis_bin}/run_orfquant.R --gtf ${gtf} --fasta ${fasta} --bam ${bam} --outdir raw \\
        --threads ${task.cpus ?: 1} --read-lengths '${params.read_lengths_orfquant}' \\
        --psite-offsets '${params.psite_offsets_orfquant}' --psite-offsets-file ${offsets}
    test -s raw/orfquant_Detected_ORFs.gtf || { echo 'ORFquant produced no output' >&2; exit 1; }
    printf '"%s":\n    ORFquant: 1.1.0\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf '##gff-version 3\nchr1\tORFquant\tCDS\t101\t220\t.\t+\t0\tORF_id=ORFQ1;transcript_id=TX1;start_codon=ATG\n' > raw/orfquant_Detected_ORFs.gtf
    printf '"stub":\n    ORFquant: stub\n' > versions.yml
    """
}
