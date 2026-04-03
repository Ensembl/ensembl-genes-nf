#!/usr/bin/env python3
"""
write_manifest.py — Write output_manifest.json for the long_read pipeline.
Required by HiveRunNextflow for dataflow back to eHive channel 2.
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pipeline', default='long_read')
    parser.add_argument('--version',  default='1.0.0')
    parser.add_argument('--outdir',   required=True)
    parser.add_argument('--gff3',     nargs='+', default=[])
    parser.add_argument('--bam',      nargs='+', default=[])
    args = parser.parse_args()

    outputs = []
    for path in args.gff3:
        p = Path(path)
        if p.exists():
            # Extract sample id from filename stem (strip .classified.gff3 etc.)
            sample_id = p.stem.split('.')[0]
            outputs.append({
                'type': 'gff3',
                'path': str(p.resolve()),
                'meta': {'id': sample_id, 'biotype': 'isoseq'}
            })

    for path in args.bam:
        p = Path(path)
        if p.exists():
            sample_id = p.stem.split('.')[0]
            outputs.append({
                'type': 'bam',
                'path': str(p.resolve()),
                'meta': {'id': sample_id}
            })

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'outputs':      outputs,
    }

    out_path = Path('output_manifest.json')
    with open(out_path, 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[write_manifest] Written {len(outputs)} outputs to {out_path}", flush=True)


if __name__ == '__main__':
    main()
