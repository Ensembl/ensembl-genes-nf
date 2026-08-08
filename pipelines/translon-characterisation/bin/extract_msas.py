#!/usr/bin/env python3
"""Extract splice-aware, instance-addressable MAF alignments.

Each attributed instance selects its chromosome's verified ``.maf.bb`` index.
No cross-chromosome fallback is attempted: a missing index is an
uninformative MSA result, not evidence against the ORF.
"""
import argparse
import json
import re
import subprocess
from pathlib import Path


def safe(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', value)


def rows(path):
    with open(path) as handle:
        yield from (json.loads(line) for line in handle if line.strip())


def fasta_has_reference(path, reference):
    return any(line.startswith(f'>{reference}') for line in Path(path).read_text().splitlines())


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--instances', required=True)
    p.add_argument('--maf-dir', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--axis-output', required=True)
    p.add_argument('--identity-output', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--maf-extract', default='mafExtract')
    p.add_argument('--stitch', default='stitch_maf_alignment.py')
    p.add_argument('--species-map', required=True)
    a = p.parse_args()
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    maf_dir = Path(a.maf_dir)
    output, identities = [], []
    for instance in rows(a.instances):
        iid, stem = instance['instance_id'], safe(instance['instance_id'])
        index = maf_dir / f"{instance['chrom']}.maf.bb"
        if not index.is_file():
            output.append({'instance_id': iid, 'axis': {'msa_extract': {'state': 'uninformative', 'reason': 'no verified chromosome MAF index', 'chrom': instance['chrom']}}})
            continue
        beds = []
        for n, block in enumerate(instance.get('exon_blocks', []), 1):
            left, right = max(instance['start'], block['start']), min(instance['end'], block['end'])
            if left < right:
                bed = outdir / f'{stem}.exon{n}.bed'
                bed.write_text(f"{instance['chrom']}\t{left}\t{right}\t{stem}_{n}\t0\t{instance['strand']}\n")
                beds.append(bed)
        if not beds:
            output.append({'instance_id': iid, 'axis': {'msa_extract': {'state': 'unattributable', 'reason': 'instance has no extractable exon blocks'}}})
            continue
        pieces = outdir / f'{stem}.maf'; pieces.mkdir(exist_ok=True)
        try:
            for n, bed in enumerate(beds, 1):
                fields = bed.read_text().rstrip().split('\t')
                region = f'{fields[0]}:{fields[1]}-{fields[2]}'
                subprocess.run([a.maf_extract, str(index), f'-region={region}', str(pieces / f'exon{n}.maf')], check=True)
            maf_files = sorted(pieces.glob('*.maf'))
            alignment = outdir / f'{stem}.msa.fasta'
            subprocess.run([a.stitch, '--mafs', *(str(item) for item in maf_files), '--strand', instance['strand'], '--species-map', a.species_map, '--output', str(alignment)], check=True)
            lines = alignment.read_text().splitlines()
            if not lines or not fasta_has_reference(alignment, a.reference):
                raise ValueError('stitched MSA lacks the declared reference sequence')
            identities.append({'instance_id': iid, 'alignment': alignment.name})
            output.append({'instance_id': iid, 'axis': {'msa_extract': {'state': 'positive', 'reason': 'splice-aware alignment extracted from verified MAF/index', 'chrom': instance['chrom'], 'alignment': alignment.name}}})
        except (OSError, subprocess.CalledProcessError, ValueError) as error:
            output.append({'instance_id': iid, 'axis': {'msa_extract': {'state': 'uninformative', 'reason': 'MAF extraction did not yield an attributable reference alignment', 'error': str(error)}}})
    Path(a.axis_output).write_text(''.join(json.dumps(item, sort_keys=True) + '\n' for item in output))
    Path(a.identity_output).write_text(''.join(json.dumps(item, sort_keys=True) + '\n' for item in identities))


if __name__ == '__main__':
    main()
