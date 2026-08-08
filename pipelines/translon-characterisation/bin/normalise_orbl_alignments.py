#!/usr/bin/env python3
"""Make ORBL_tools downloads safe per-instance PhyloCSF inputs.

The official downloader does not promise a stable filename convention.  The
request file is therefore the identity source.  A downloaded file is accepted
only when exactly one request identifier can be matched; ambiguous or missing
files are retained in a machine-readable report rather than silently scored.
"""
import argparse
import json
import re
import shutil
from pathlib import Path


def safe(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', value)


def requests(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) >= 3:
            result[fields[2]] = fields[2]
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--requests', required=True)
    p.add_argument('--indir', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--report', required=True)
    a = p.parse_args()
    ids = requests(a.requests)
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    report = []
    for source in sorted(Path(a.indir).glob('*.fa*')):
        matches = [iid for iid in ids if iid in source.name or safe(iid) in source.name]
        if len(matches) != 1:
            report.append({'source': source.name, 'state': 'unattributable', 'reason': 'downloaded alignment filename does not map uniquely to an ORBL request', 'matches': matches})
            continue
        iid = matches[0]
        destination = outdir / f'{safe(iid)}.fa'
        lines = source.read_text().splitlines()
        if not lines or not lines[0].startswith('>'):
            report.append({'source': source.name, 'instance_id': iid, 'state': 'uninformative', 'reason': 'downloaded alignment is not FASTA'})
            continue
        header = lines[0][1:].split()[0]
        lines[0] = f'>{header} instance_id={iid}'
        destination.write_text('\n'.join(lines) + '\n')
        report.append({'source': source.name, 'instance_id': iid, 'output': destination.name, 'state': 'positive'})
    Path(a.report).write_text(''.join(json.dumps(row, sort_keys=True) + '\n' for row in report))


if __name__ == '__main__':
    main()
