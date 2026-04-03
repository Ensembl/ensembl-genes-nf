#!/usr/bin/env python3
"""write_manifest.py — Write output_manifest.json for HiveRunNextflow."""

from __future__ import annotations

import argparse, json
from datetime import datetime, timezone
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pipeline', required=True)
    parser.add_argument('--version',  required=True)
    parser.add_argument('--outdir',   required=True)
    parser.add_argument('--gff3',     required=True, nargs='+')
    args = parser.parse_args()

    outputs = []
    for gff3 in args.gff3:
        p = Path(gff3)
        if p.exists():
            outputs.append({'type': 'ncrna_gff3', 'path': str(p.resolve()),
                            'meta': {'size_bytes': p.stat().st_size}})

    manifest = {
        'pipeline': args.pipeline, 'version': args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs': outputs,
    }
    out_path = Path(args.outdir) / 'output_manifest.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(f'[write_manifest] Wrote {out_path}', flush=True)


if __name__ == '__main__':
    main()
