process CANONICALISE_PACBIO_READS {
    tag "${meta.id}:${meta.classification}"
    label 'process_high_memory'
    container { params.ccs_container ?: 'https://depot.galaxyproject.org/singularity/pbccs:6.0.0--h9ee0642_0' }

    input:
    tuple val(meta), path(input_bam)
    path inspector

    output:
    tuple val(meta), path('canonical.fastq.gz'), emit: reads
    tuple val(meta), path('molecule_audit.tsv'), emit: molecule_audit
    path 'versions.yml', emit: versions

    script:
    def ccs_args = task.ext.ccs_args ?: params.ccs_args ?: ''
    """
    set -euo pipefail
    input_bam=\$(find . -maxdepth 1 -name '*.bam' -print -quit)
    input_records=\$(samtools view -c -F 0x900 "\${input_bam}")
    if [ "${meta.classification}" = "PACBIO_SUBREAD_BAM" ]; then
      ${params.ccs_command ?: 'ccs'} ${ccs_args} "\${input_bam}" consensus.bam
      ccs_version=\$(${params.ccs_command ?: 'ccs'} --version 2>&1 | head -n1)
      test -n "${params.ccs_expected_version ?: ''}"
      echo "\${ccs_version}" | grep -F "${params.ccs_expected_version ?: ''}" >/dev/null
      samtools fastq -F 0x900 consensus.bam | gzip -c > canonical.fastq.gz
    else
      ccs_version='not_run'
      samtools fastq -F 0x900 "\${input_bam}" | gzip -c > canonical.fastq.gz
    fi
    python3 ${inspector} validate-fastq canonical.fastq.gz PACBIO_CCS validation.json
    canonical_sha256=\$(sha256sum canonical.fastq.gz | awk '{print \$1}')
    printf 'run_accession\\tclassification\\tinput_records\\toutput_records\\tdistinct_molecules\\tcanonical_fastq_sha256\\tccs_version\\n' > molecule_audit.tsv
    printf '${meta.id}\\t${meta.classification}\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
      "\${input_records}" \\
      "\$(sed -n 's/.*\"records\": \([0-9]*\).*/\1/p' validation.json)" \\
      "\$(sed -n 's/.*\"distinct_molecules\": \([0-9]*\).*/\1/p' validation.json)" \\
      "\${canonical_sha256}" \\
      "\${ccs_version}" >> molecule_audit.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ccs: \${ccs_version}
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    printf '@movie/1/ccs\\nACGT\\n+\\n!!!!\\n' | gzip -c > canonical.fastq.gz
    printf 'run_accession\\tclassification\\tinput_records\\toutput_records\\tdistinct_molecules\\tcanonical_fastq_sha256\\tccs_version\\n${meta.id}\\t${meta.classification}\\t1\\t1\\t1\\tstub\\tstub\\n' > molecule_audit.tsv
    printf '"%s":\\n    ccs: stub\\n    samtools: 1.20\\n' '${task.process}' > versions.yml
    """
}
