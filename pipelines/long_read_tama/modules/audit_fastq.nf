process AUDIT_FASTQ {
    tag "${meta.id}:${meta.classification}"
    label 'process_light'

    input:
    tuple val(meta), path(probe), path(validation)

    output:
    tuple val(meta), path('molecule_audit.tsv'), emit: molecule_audit
    path 'versions.yml', emit: versions

    script:
    """
    python3 - <<'PY'
import json
from pathlib import Path

probe = json.loads(Path('${probe}').read_text())
validation = json.loads(Path('${validation}').read_text())
fields = ['run_accession', 'classification', 'header_representation', 'records_sampled',
          'distinct_ids', 'distinct_molecules', 'malformed', 'status']
values = [
    '${meta.id}', '${meta.classification ?: 'UNCLASSIFIED'}', probe.get('header_representation', 'UNKNOWN'),
    probe.get('records_sampled', 0), validation.get('distinct_ids', 0),
    validation.get('distinct_molecules', 0), probe.get('malformed_count', 0),
    validation.get('status', 'VALIDATED'),
]
Path('molecule_audit.tsv').write_text('\\t'.join(fields) + '\\n' + '\\t'.join(map(str, values)) + '\\n')
PY
    printf '"%s":\n    audit: generated\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'run_accession\tclassification\theader_representation\trecords_sampled\tdistinct_ids\tdistinct_molecules\tmalformed\tstatus\n${meta.id}\t${meta.classification ?: 'UNCLASSIFIED'}\tUNKNOWN\t1\t1\t1\t0\tstub\n' > molecule_audit.tsv
    printf '"%s":\n    audit: generated\n' '${task.process}' > versions.yml
    """
}
