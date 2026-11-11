process ACQUIRE_READS {
    tag "${meta.sample_id}"
    container 'quay.io/biocontainers/fastq-dl:4.0.1--pyhdfd78af_0'
    label 'cage_python'
    input:
    tuple val(meta), val(record)

    output:
    tuple val(meta), path('*.fastq*'), emit: reads

    script:
    def source = record.fastq_1 ? "cp -L ${record.fastq_1} ${record.sample_id}_R1.fastq${record.fastq_2 ? " && cp -L ${record.fastq_2} ${record.sample_id}_R2.fastq" : ''}" : "fastq-dl --accession ${record.accession} --outdir ."
    def checksum_1 = record.fastq_1 && record.checksum_1 && record.checksum_1 != '-' ? "echo '${record.checksum_1}  ${record.sample_id}_R1.fastq' | sha256sum -c -" : ':'
    def checksum_2 = record.fastq_2 && record.checksum_2 && record.checksum_2 != '-' ? "echo '${record.checksum_2}  ${record.sample_id}_R2.fastq' | sha256sum -c -" : ':'
    """
    ${source}
    if [[ ! -s ${record.sample_id}_R1.fastq && -s ${record.accession ?: 'missing.fastq'} ]]; then cp -L ${record.accession} ${record.sample_id}_R1.fastq; fi
    test -n "\$(find . -maxdepth 1 -type f -name '*.fastq*' -print -quit)" || { echo 'No FASTQ files acquired' >&2; exit 1; }
    ${checksum_1}
    ${checksum_2}
    """

    stub:
    """
    printf '@${meta.sample_id}\\nACGT\\n+\\nIIII\\n' > ${meta.sample_id}_R1.fastq
    """
}
