/* The stable caller contract is RUN_<TOOL> -> STANDARDISE_<TOOL>.
 * Preparation is deliberately inside the runner so each caller has only two
 * externally visible processes. */

process RUN_RIBOCODE {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda 'bioconda::ribocode'
    container 'quay.io/biocontainers/ribocode:1.2.15--pyhdc42f0e_1'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_ribocode ?: ''
    def threads = task.cpus ?: 1
    """
    mkdir -p raw prep
    export MPLCONFIGDIR=\$PWD/.mplconfig
    mkdir -p "\$MPLCONFIGDIR"
    python ${projectDir}/bin/make_transcriptome_annotation.py \\
        --gtf ${gtf} --fasta ${fasta} \\
        --out-gtf raw/transcriptome.gtf --out-fasta raw/transcriptome.fa
    RiboCode_onestep -g raw/transcriptome.gtf -f raw/transcriptome.fa -r ${bam} -l no -s ATG -o raw/ribocode -t ${threads} ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'RiboCode produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RiboCode: \$(RiboCode_onestep --version 2>&1 | head -n1 || true)
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'chrom\tstart\tend\ttranscript_id\tframe\tscore\tpval\nchr1\t100\t200\tTX1\t0\t10\t0.01\n' > raw/ribocode.tsv
    printf '"stub":\n    RiboCode: stub\n' > versions.yml
    """
}

process RUN_RIBOTRICER {
    tag "${meta.id}"
    label 'process_single_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda 'bioconda::ribotricer'
    container 'quay.io/biocontainers/ribotricer:1.5.0--pyhdfd78af_0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_ribotricer ?: ''
    """
    mkdir -p raw prep
    export MPLCONFIGDIR=\$PWD/.mplconfig
    mkdir -p "\$MPLCONFIGDIR"
    ribotricer prepare-orfs --gtf ${gtf} --fasta ${fasta} --prefix raw/orfs --start_codons ATG
    ribotricer detect-orfs --bam ${bam} --ribotricer_index raw/orfs_candidate_orfs.tsv --prefix raw/ribotricer ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'Ribotricer produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ribotricer: \$(ribotricer --version 2>&1 | head -n1 || true)
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'transcript_id\tstart\tend\tframe\tscore\tperiodicity\nTX1\t120\t210\t0\t5.2\t0.8\n' > raw/ribotricer.tsv
    printf '"stub":\n    Ribotricer: stub\n' > versions.yml
    """
}

process RUN_ORFQUANT {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfquant:1.1.0-txdbmaker2'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai), path(offsets)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    """
    export MPLCONFIGDIR=\$PWD/.mplconfig
    mkdir -p raw "\$MPLCONFIGDIR"
    Rscript ${projectDir}/bin/run_orfquant.R \\
        --gtf ${gtf} --fasta ${fasta} --bam ${bam} --outdir raw \\
        --threads ${task.cpus ?: 1} \\
        --read-lengths '${params.read_lengths_orfquant}' \\
        --psite-offsets '${params.psite_offsets_orfquant}' \\
        --psite-offsets-file ${offsets}
    test -s raw/orfquant_Detected_ORFs.gtf || { echo 'ORFquant produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ORFquant: \$(Rscript -e 'cat(as.character(packageVersion("ORFquant")))')
        R: \$(Rscript --version 2>&1 | head -n1)
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'transcript_id\torf_id\tstart\tend\treads\tTPM\nTX1\tORFQ1\t100\t220\t120\t3.2\n' > raw/orfquant.tsv
    printf '"stub":\n    ORFquant: stub\n' > versions.yml
    """
}

process RUN_RPBP {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(fastq)
    path gtf
    path fasta
    path ribosomal_fasta
    path adapter_fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def adapter_yaml = "adapter_file: ${adapter_fasta}"
    """
    mkdir -p raw
    cat > raw/rpbp.yaml <<-END_CONFIG
    gtf: ${gtf}
    fasta: ${fasta}
    ribosomal_fasta: ${ribosomal_fasta}
    genome_name: ${meta.id}_rpbp
    genome_base_path: raw/genome
    ribosomal_index: raw/ribosomal_index
    star_index: raw/star_index
    riboseq_samples:
      ${meta.id}: ${fastq}
    ${adapter_yaml}
    riboseq_data: raw/data
    END_CONFIG
    prepare-rpbp-genome raw/rpbp.yaml --star-options '--genomeSAindexNbases 8' --num-cpus ${task.cpus ?: 1} --logging-level INFO
    run-all-rpbp-instances raw/rpbp.yaml --profiles-only --num-cpus ${task.cpus ?: 1} --logging-level INFO
    test -n "\$(find raw -type f | head -1)" || { echo 'Rp-Bp produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rp-Bp: 3.0.1
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'chr1\t300\t420\tORF_RP1\t+\t8.5\n' > raw/rpbp.bed
    printf '"stub":\n    Rp-Bp: stub\n' > versions.yml
    """
}

process RUN_IRIBO {
    tag "${meta.id}"
    label 'process_ultra_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-iribo:1.0.0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_iribo ?: ''
    def threads = task.cpus ?: 2
    """
    mkdir -p raw
    printf '%s\n' '${bam}' > raw/riboseq_bams.txt
    iRibo --RunMode=GetCandidateORFs --Genome=${fasta} --Annotations=${gtf} --Output=raw/candidates --Threads=${threads} ${args}
    iRibo --RunMode=GenerateTranslationProfile --Genome=${fasta} --Annotations=${gtf} --Riboseq=raw/riboseq_bams.txt --CandidateORFs=raw/candidates/candidate_orfs --Output=raw/profile --Threads=${threads} ${args}
    # GenerateTranslationProfile defaults to one scramble; pass the same
    # value here because GenerateTranslatome.R otherwise defaults to 100 and
    # indexes non-existent scrambled1..scrambled100 columns.
    Rscript \$IRIBO_HOME/GenerateTranslatome.R --TranslationCalls=raw/profile/translation_calls --NullDistribution=raw/profile/null_distribution --CandidateORFs=raw/candidates/candidate_orfs --Output=raw/translatome --Threads=${threads} --Scrambles=1
    test -s raw/translatome/translated_orfs.csv || { echo 'iRibo produced no native translated_orfs.csv' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iRibo: source-pinned
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'stub\n' > raw/iribo_stub.tsv
    printf '"stub":\n    iRibo: stub\n' > versions.yml
    """
}

process RUN_ORFRATER {
    tag "${meta.id}"
    label 'process_single_high_memory'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfrater:1.0.0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai), path(bed12)
    path gtf
    path fasta
    path model
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_orfrater ?: ''
    """
    mkdir -p raw
    python \$ORFRATER_HOME/find_orfs_and_types.py --help >/dev/null
    test -n "\$(find ${model} -type f | head -1)" || { echo 'ORF-RATER model bundle is empty' >&2; exit 2; }
    python \$ORFRATER_HOME/quantify_orfs.py ${bam} --inbed ${bed12} --subdir raw --ratingsfile ${model}/orfratings.h5 --metagenefile ${model}/metagene.txt --offsetfile ${model}/offsets.txt ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'ORF-RATER produced no native output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ORF-RATER: source-pinned
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'stub\n' > raw/orfrater_stub.tsv
    printf '"stub":\n    ORF-RATER: stub\n' > versions.yml
    """
}

process RUN_RIBORF {
    tag "${meta.id}"
    label 'process_single_high_memory'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-riborf:1.0.0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(reads_sam), path(genepred), path(bed12)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_riborf ?: ''
    """
    mkdir -p raw
    perl \$RIBORF_HOME/ribORF.pl -f ${reads_sam} -c ${genepred} -o raw/predictions ${args}
    test -n "\$(find raw/predictions -type f | head -1)" || { echo 'RibORF produced no native output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RibORF: 2.0
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'stub\n' > raw/riborf_stub.tsv
    printf '"stub":\n    RibORF: stub\n' > versions.yml
    """
}

process RUN_RIBOTISH {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/ribotish:0.2.8--pyhdfd78af_0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_ribotish ?: ''
    def threads = task.cpus ?: 1
    """
    mkdir -p raw
    ribotish predict -b ${bam} -g ${gtf} -f ${fasta} -o raw/ribotish.txt --blocks --seq --aaseq -p ${threads} ${args}
    test -s raw/ribotish.txt || { echo 'Ribo-TISH produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RiboTISH: \$(ribotish --version 2>&1 | head -n1 || true)
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'stub\n' > raw/ribotish_stub.tsv
    printf '"stub":\n    RiboTISH: stub\n' > versions.yml
    """
}

process RUN_RIBOTIE {
    tag "${meta.id}"
    label 'process_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container params.ribotie_gpu ? \
        'ghcr.io/jackcurragh/translon-ribotie-cuda:1.0.0' : \
        'ghcr.io/jackcurragh/translon-ribotie:1.0.0'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'raw', saveAs: { filename -> "${meta.id}/${task.process}/raw" }
    publishDir "${params.outdir}/native_outputs", mode: 'copy', pattern: 'versions.yml', saveAs: { filename -> "${meta.id}/${task.process}/${filename}" }
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: params.args_ribotie ?: ''
    """
    mkdir -p raw
    python - <<'PY'
    import torch
    if not torch.cuda.is_available():
        raise SystemExit('RiboTIE requires CUDA, but torch.cuda.is_available() is false')
    print('RiboTIE CUDA device:', torch.cuda.get_device_name(0))
    PY
    cat > raw/ribotie.yml <<-END_CONFIG
    gtf_path: ${gtf}
    fa_path: ${fasta}
    ribo_paths:
      ${meta.id}: ${bam}
    h5_path: raw/${meta.id}.h5
    END_CONFIG
    ribotie raw/ribotie.yml ${args}
    test -n "\$(find raw -type f \( -name '*.csv' -o -name '*.gtf' \) | head -1)" || { echo 'RiboTIE produced no native result table' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RiboTIE: \$(ribotie --version 2>&1 | head -n1 || true)
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'stub\n' > raw/ribotie_stub.tsv
    printf '"stub":\n    RiboTIE: stub\n' > versions.yml
    """
}

process STANDARDISE_CALLER {
    tag "${meta.id}"
    label 'process_light'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda 'conda-forge::python=3.11'
    // The standardiser is Python-only, but Nextflow task metrics also need
    // `ps` inside the container; the pinned tool image provides both.
    container 'quay.io/biocontainers/ribotricer:1.5.0--pyhdfd78af_0'
    input:
    tuple val(meta), path(rawdir)
    val tool
    path gtf
    output:
    tuple val(meta), val(tool), path("${meta.id}_${tool}.bed12"), emit: bed12
    tuple val(meta), val(tool), path("${meta.id}_${tool}.tsv"), emit: standardized
    path 'versions.yml', emit: versions, topic: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/bin/standardise_caller.py --adapter-version 2 --raw ${rawdir} --tool ${tool} --sample ${meta.id} --gtf ${gtf} --output ${meta.id}_${tool}.tsv --bed12 ${meta.id}_${tool}.bed12 ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardiser: python3
    END_VERSIONS
    """
    stub:
    """
    printf 'sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\ttranscript_id\torf_id\tscore\tpval\tqval\textra_json\n${meta.id}\t${tool}\tchr1\t100\t200\t+\t0\tTX1\tORF1\t10\t0.01\t\t{}\n' > ${meta.id}_${tool}.tsv
    printf 'chr1\t100\t200\tORF1|TX1|${tool}\t900\t+\t100\t200\t0,0,0\t1\t100,\t0,\n' > ${meta.id}_${tool}.bed12
    printf '"stub":\n    standardiser: stub\n' > versions.yml
    """
}
