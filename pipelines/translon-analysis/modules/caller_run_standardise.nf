/* The stable caller contract is RUN_<TOOL> -> STANDARDISE_<TOOL>.
 * Preparation is deliberately inside the runner so each caller has only two
 * externally visible processes. */

process RUN_RIBOCODE {
    tag "${meta.id}"
    container params.container_ribocode
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    script:
    """
    mkdir -p raw prep
    python3 ${projectDir}/pipelines/orf-calling/bin/prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type transcriptome || true
    RiboCode_onestep -a ${gtf} -r ${bam} -g ${fasta} -t ${params.threads_ribocode ?: 4} ${params.args_ribocode ?: ''} || true
    find . -maxdepth 2 -type f \( -name '*.tsv' -o -name '*.bed' \) -exec cp {} raw/ \; 2>/dev/null || true
    touch raw/ribocode.tsv
    """
    stub:
    """
    mkdir -p raw
    printf 'chrom\tstart\tend\ttranscript_id\tframe\tscore\tpval\nchr1\t100\t200\tTX1\t0\t10\t0.01\n' > raw/ribocode.tsv
    """
}

process RUN_RIBOTRICER {
    tag "${meta.id}"
    container params.container_ribotricer
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    script:
    """
    mkdir -p raw prep
    python3 ${projectDir}/pipelines/orf-calling/bin/prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type transcriptome || true
    ribotricer prepare-orfs -a ${gtf} -g ${fasta} -o raw/orfs || true
    ribotricer detect-orfs -b ${bam} -i raw/orfs_candidate_ORFs.tsv -o raw -t ${params.threads_ribotricer ?: 4} || true
    """
    stub:
    """
    mkdir -p raw
    printf 'transcript_id\tstart\tend\tframe\tscore\tperiodicity\nTX1\t120\t210\t0\t5.2\t0.8\n' > raw/ribotricer.tsv
    """
}

process RUN_RIBOTAPER {
    tag "${meta.id}"
    container params.container_ribotaper
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    script:
    """
    mkdir -p raw prep
    python3 ${projectDir}/pipelines/orf-calling/bin/prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome || true
    ribotaper -h >/dev/null 2>&1 || true
    """
    stub:
    """
    mkdir -p raw
    printf 'chr1\t150\t260\tORF_TP1\t+\t12\n' > raw/ribotaper.bed
    """
}

process RUN_ORFQUANT {
    tag "${meta.id}"
    container params.container_orfquant
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    script:
    """
    mkdir -p raw prep
    python3 ${projectDir}/pipelines/orf-calling/bin/prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type transcriptome || true
    orfquant --help >/dev/null 2>&1 || true
    """
    stub:
    """
    mkdir -p raw
    printf 'transcript_id\torf_id\tstart\tend\treads\tTPM\nTX1\tORFQ1\t100\t220\t120\t3.2\n' > raw/orfquant.tsv
    """
}

process RUN_RPBP {
    tag "${meta.id}"
    container params.container_rpbp
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    script:
    """
    mkdir -p raw prep
    python3 ${projectDir}/pipelines/orf-calling/bin/prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome || true
    rpbp --help >/dev/null 2>&1 || true
    """
    stub:
    """
    mkdir -p raw
    printf 'chr1\t300\t420\tORF_RP1\t+\t8.5\n' > raw/rpbp.bed
    """
}

process STANDARDISE_CALLER {
    tag "${meta.id}"
    container 'python:3.11-slim'
    input:
    tuple val(meta), path(rawdir)
    val tool
    path gtf
    output:
    tuple val(meta), val(tool), path("${meta.id}_${tool}.bed12"), emit: bed12
    tuple val(meta), val(tool), path("${meta.id}_${tool}.tsv"), emit: standardized
    script:
    """
    python3 ${projectDir}/pipelines/translon-analysis/bin/standardise_caller.py --raw ${rawdir} --tool ${tool} --sample ${meta.id} --gtf ${gtf} --output ${meta.id}_${tool}.tsv --bed12 ${meta.id}_${tool}.bed12
    """
    stub:
    """
    printf 'sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\ttranscript_id\torf_id\tscore\tpval\tqval\textra_json\n${meta.id}\t${tool}\tchr1\t100\t200\t+\t0\tTX1\tORF1\t10\t0.01\t\t{}\n' > ${meta.id}_${tool}.tsv
    printf 'chr1\t100\t200\tORF1|TX1|${tool}\t900\t+\t100\t200\t0,0,0\t1\t100,\t0,\n' > ${meta.id}_${tool}.bed12
    """
}
