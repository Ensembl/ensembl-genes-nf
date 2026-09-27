#!/usr/bin/env python3
"""Run set-level product clustering without turning absent structures into data.

MMseqs2 clusters trusted translated peptide products.  Foldseek is attempted
only for records with explicitly supplied, readable structure paths; otherwise
the structural dimension is recorded as uninformative for each product.
"""
import argparse
import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def rows(path):
    with open(path) as handle:
        yield from (json.loads(line) for line in handle if line.strip())


def clusters_from_tsv(path):
    result = defaultdict(list)
    if not path.exists():
        return result
    for line in path.read_text().splitlines():
        representative, member = line.split('\t')[:2]
        result[member].append(representative)
    return result


def run_cluster(executable, fasta, output_prefix, workdir, extra=()):
    subprocess.run([executable, 'easy-cluster', str(fasta), str(output_prefix), str(workdir), *extra], check=True)
    return clusters_from_tsv(Path(f'{output_prefix}_cluster.tsv'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--instances', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--mmseqs', default='mmseqs')
    p.add_argument('--foldseek', default='foldseek')
    p.add_argument('--min-seq-id', default='0.9')
    a = p.parse_args()
    products = {}
    for item in rows(a.instances):
        sequence = item.get('peptide_sequence')
        if sequence:
            products.setdefault(item.get('peptide_id') or f"peptide:{sequence}", {'sequence': sequence, 'instances': []})['instances'].append(item['instance_id'])
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fasta = root / 'products.fa'
        fasta.write_text(''.join(f'>{pid}\n{value["sequence"]}\n' for pid, value in products.items()))
        sequence_clusters = run_cluster(a.mmseqs, fasta, root / 'mmseqs', root / 'mmseqs_tmp', ('--min-seq-id', str(a.min_seq_id), '-c', '0.8')) if products else {}
        structures = {pid: Path(path) for pid, value in products.items() if (path := next((x.get('structure_path') for x in rows(a.instances) if (x.get('peptide_id') or f"peptide:{x.get('peptide_sequence', '')}") == pid and x.get('structure_path')), None)) and Path(path).is_file()}
        structural_clusters, structural_error = {}, None
        if structures:
            structure_list = root / 'structures.txt'; structure_list.write_text('\n'.join(str(path) for path in structures.values()) + '\n')
            try:
                structural_clusters = run_cluster(a.foldseek, structure_list, root / 'foldseek', root / 'foldseek_tmp')
            except (OSError, subprocess.CalledProcessError) as error:
                structural_error = str(error)
    with open(a.output, 'w') as out:
        for pid, value in products.items():
            record = {'peptide_id': pid, 'instance_ids': value['instances'], 'sequence_clusters': sequence_clusters.get(pid, []), 'sequence_cluster_state': 'positive' if sequence_clusters.get(pid) else 'uninformative'}
            if pid not in structures:
                record['structural_cluster'] = {'state': 'uninformative', 'reason': 'no attributable structure model supplied'}
            elif structural_error:
                record['structural_cluster'] = {'state': 'uninformative', 'reason': 'Foldseek did not complete', 'error': structural_error}
            else:
                record['structural_cluster'] = {'state': 'positive' if structural_clusters.get(pid) else 'null', 'representatives': structural_clusters.get(pid, [])}
            out.write(json.dumps(record, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
