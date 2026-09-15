process PREPARE_LONG_READ_INPUT {
    tag "${meta.id}"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(meta), path(reads)
    path inspector

    output:
    tuple val(meta), path('reads.fastq.gz'), emit: reads
    tuple val(meta), path('read_validation.tsv'), emit: report
    tuple val(meta), path('molecule_audit.tsv'), emit: molecule_audit
    path 'versions.yml', emit: versions

    script:
    """
    test \$(find ${reads} -maxdepth 0 -type f | wc -l) -eq 1
    gzip -t ${reads}
    seqkit stats -T ${reads} > read_validation.tsv
    python3 ${inspector} probe ${reads} header_probe.json
    if [ "${meta.classification}" != "UNCLASSIFIED" ]; then
      python3 ${inspector} validate-fastq ${reads} ${meta.expected_header_representation} full_validation.json
    fi
    printf 'run\\tclassification\\theader_representation\\trecords_sampled\\tmalformed\\tstatus\\n' > molecule_audit.tsv
    printf '${meta.id}\\t${meta.classification ?: "UNCLASSIFIED"}\\t%s\\t%s\\t%s\\t%s\\n' \\
      "\$(sed -n 's/.*\"header_representation\": \"\\([^\"]*\\)\".*/\\1/p' header_probe.json)" \\
      "\$(sed -n 's/.*\"records_sampled\": \\([0-9]*\\).*/\\1/p' header_probe.json)" \\
      "\$(sed -n 's/.*\"malformed_count\": \\([0-9]*\\).*/\\1/p' header_probe.json)" \\
      "\$(test -s header_probe.json && echo ok || echo failed)" >> molecule_audit.tsv
    printf '${meta.id}\\t${meta.classification ?: "UNCLASSIFIED"}\\t%s\\t%s\\t%s\\t%s\\n' \\
      "\$(sed -n 's/.*\"header_representation\": \"\\([^\"]*\\)\".*/\\1/p' header_probe.json)" \\
      "\$(sed -n 's/.*\"records_sampled\": \\([0-9]*\\).*/\\1/p' header_probe.json)" \\
      "\$(sed -n 's/.*\"malformed_count\": \\([0-9]*\\).*/\\1/p' header_probe.json)" \\
      "\$(test -s header_probe.json && echo ok || echo failed)" >> molecule_audit.tsv
    if [ "${meta.classification}" != "UNCLASSIFIED" ]; then
      printf '%s\\n' "\$(sed -n 's/.*\"records\": \\([0-9]*\\).*/\\1/p' full_validation.json)" > /dev/null
    fi
    cp ${reads} reads.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    printf '@stub\\nACGT\\n+\\n!!!!\\n' | gzip -c > reads.fastq.gz
    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\nreads.fastq.gz\\tFASTQ\\tDNA\\t1\\t4\\t4\\t4\\t4\\n' > read_validation.tsv
    printf 'run\\tclassification\\theader_representation\\trecords_sampled\\tdistinct_ids\\tdistinct_molecules\\tmalformed\\tstatus\\n${meta.id}\\t${meta.classification ?: "UNCLASSIFIED"}\\tUNKNOWN\\t1\\t1\\t1\\t0\\tstub\\n' > molecule_audit.tsv
    printf '"%s":\\n    seqkit: 2.8.2\\n' '${task.process}' > versions.yml
    """
}
