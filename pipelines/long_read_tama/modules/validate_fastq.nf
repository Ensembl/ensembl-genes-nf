process VALIDATE_FASTQ {
    tag "${meta.id}:${meta.classification}"
    label 'process_light'
    // The staged inspector and the FASTQ statistics helper both require Python.
    // Keep this process on a runtime that actually provides python3.
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(reads)
    path inspector

    output:
    tuple val(meta), path('reads.fastq.gz'), emit: reads
    tuple val(meta), path('read_validation.tsv'), emit: report
    tuple val(meta), path('molecule_audit.tsv'), emit: molecule_audit
    path 'versions.yml', emit: versions

    script:
    def classification = meta.classification ?: 'UNCLASSIFIED'
    def expected = meta.expected_header_representation ?: 'UNKNOWN'
    """
    test -f "${reads}" || { echo "Expected FASTQ file for ${meta.id}" >&2; exit 1; }
    gzip -t "${reads}" || { echo "Invalid gzip FASTQ for ${meta.id}: ${reads}" >&2; exit 1; }
    python3 "${inspector}" stats "${reads}" read_validation.tsv || { echo "FASTQ structure validation failed for ${meta.id}" >&2; exit 1; }
    python3 "${inspector}" probe "${reads}" header_probe.json
    if [ "${classification}" != "UNCLASSIFIED" ]; then
        python3 "${inspector}" validate-fastq "${reads}" "${expected}" full_validation.json
    fi
    test -s header_probe.json || { echo "Header probe produced no result for ${meta.id}" >&2; exit 1; }
    python3 - <<'PY'
import json
from pathlib import Path

probe = json.loads(Path('header_probe.json').read_text())
validation = json.loads(Path('full_validation.json').read_text()) if Path('full_validation.json').exists() else {}
fields = ['run_accession', 'classification', 'header_representation', 'records_sampled',
          'distinct_ids', 'distinct_molecules', 'malformed', 'status']
values = [
    '${meta.id}', '${classification}', probe.get('header_representation', 'UNKNOWN'),
    probe.get('records_sampled', 0), validation.get('distinct_ids', probe.get('distinct_ids', 0)),
    validation.get('distinct_molecules', probe.get('distinct_molecules', 0)),
    validation.get('malformed_count', probe.get('malformed_count', 0)),
    validation.get('status', 'PROBE_ONLY'),
]
Path('molecule_audit.tsv').write_text('\t'.join(fields) + '\n' + '\t'.join(map(str, values)) + '\n')
PY
    cp -p "${reads}" reads.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    printf '@stub/ccs\nACGT\n+\n!!!!\n' | gzip -c > reads.fastq.gz
    printf 'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\nreads.fastq.gz\tFASTQ\tDNA\t1\t4\t4\t4\t4\n' > read_validation.tsv
    printf 'run_accession\tclassification\theader_representation\trecords_sampled\tdistinct_ids\tdistinct_molecules\tmalformed\tstatus\n${meta.id}\t${meta.classification ?: 'UNCLASSIFIED'}\tUNKNOWN\t1\t1\t1\t0\tstub\n' > molecule_audit.tsv
    printf '"%s":\n    python: 3.11.0\n' '${task.process}' > versions.yml
    """
}
