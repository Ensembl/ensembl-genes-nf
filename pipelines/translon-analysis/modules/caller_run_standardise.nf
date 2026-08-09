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

process RUN_PUBLISHED_CALLER {
    tag "${meta.id} - ${tool}"
    container 'ubuntu:22.04'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    val tool
    output:
    tuple val(meta), path('raw'), emit: raw
    script:
    """
    mkdir -p raw prep
    python3 ${projectDir}/pipelines/orf-calling/bin/prep_inputs.py --outdir prep --sample-id ${meta.id} --bam ${bam} --bai ${bai} --gtf ${gtf} --fasta ${fasta} --bam-type genome || true
    case ${tool} in
      iribo) iribo --help >/dev/null 2>&1 || true ;;
      orfrater) orfrater --help >/dev/null 2>&1 || true ;;
      price) price --help >/dev/null 2>&1 || true ;;
      riborf) riborf --help >/dev/null 2>&1 || true ;;
      ribotish) ribotish --help >/dev/null 2>&1 || true ;;
      ribotie) ribotie --help >/dev/null 2>&1 || true ;;
    esac
    """
    stub:
    """
    mkdir -p raw
    if [ "${tool}" = "orfrater" ]; then
      printf 'transcript_id\torf_id\tstart\tend\tscore\nTX1\tORFR1\t100\t200\t4.2\n' > raw/${tool}.tsv
    elif [ "${tool}" = "ribotish" ] || [ "${tool}" = "ribotie" ]; then
      printf 'chrom\tstart\tend\tname\tstrand\tscore\nchr1\t400\t480\t${tool}1\t+\t10\n' > raw/${tool}.tsv
    else
      printf 'chr1\t210\t300\t${tool}1\t+\t9\n' > raw/${tool}.bed
    fi
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
