#!/usr/bin/env python3
"""
write_manifest.py — Write output_manifest.json for the HiveRunNextflow bridge.

Usage:
    write_manifest.py --pipeline igtr --version 1.0.0
                      --outdir /path/to/outdir --gff3 igtr.gff3
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline', required=True)
    parser.add_argument('--version',  required=True)
    parser.add_argument('--outdir',   required=True)
    parser.add_argument('--gff3',     required=True, nargs='+',
                        help='Clustered GFF3 output file(s)')
    args = parser.parse_args()

    outputs = []
    for gff3 in args.gff3:
        p = Path(gff3)
        if p.exists():
            outputs.append({
                'type': 'igtr_gff3',
                'path': str(p.resolve()),
                'meta': {'size_bytes': p.stat().st_size},
            })

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs':      outputs,
    }

    out_path = Path(args.outdir) / 'output_manifest.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(manifest, fh, indent=2)

    print(f'[write_manifest] Wrote {out_path} ({len(outputs)} outputs)', flush=True)


if __name__ == '__main__':
    main()
