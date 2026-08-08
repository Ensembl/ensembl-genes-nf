#!/usr/bin/env python3
"""Deduplicate translated products by (interval_id, resolved frame).

The translation producer must provide `peptide_sequence`; if it cannot, the
record is retained with an unattributable peptide axis rather than translated
from a guessed coordinate/frame.
"""
import argparse, json

def read(path):
    with open(path) as h:
        yield from (json.loads(line) for line in h if line.strip())

def main():
    p = argparse.ArgumentParser(); p.add_argument('--instances', required=True); p.add_argument('--peptides', required=True); p.add_argument('--fanback', required=True); a = p.parse_args()
    products, links = {}, []
    for instance in read(a.instances):
        frame = instance.get('translation', {}).get('frame', instance.get('frame', '.'))
        key = f"{instance['interval_id']}|{frame}"
        sequence = instance.get('peptide_sequence', '')
        product = products.setdefault(key, {'peptide_id': key, 'interval_id': instance['interval_id'], 'frame': frame, 'sequence': sequence})
        if not sequence:
            product['peptide_admissibility'] = {'state': 'unattributable', 'reason': 'no trusted translated peptide sequence supplied'}
        links.append({'instance_id': instance['instance_id'], 'peptide_id': key})
    with open(a.peptides, 'w') as h:
        for row in products.values(): h.write(json.dumps(row, sort_keys=True) + '\n')
    with open(a.fanback, 'w') as h:
        for row in links: h.write(json.dumps(row, sort_keys=True) + '\n')
if __name__ == '__main__': main()
