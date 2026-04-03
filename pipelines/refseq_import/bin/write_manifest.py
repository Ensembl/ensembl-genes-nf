#!/usr/bin/env python3
"""Write output_manifest.json for HiveRunNextflow dataflow."""

import argparse
import json
import os
from datetime import datetime, timezone


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pipeline', default='refseq_import')
    ap.add_argument('--version',  default='1.0.0')
    ap.add_argument('--outdir',   required=True)
    ap.add_argument('--gff3',     required=True)
    args = ap.parse_args()

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs': [
            {
                'type': 'gff3',
                'path': os.path.join(args.outdir, 'refseq_import', args.gff3),
            }
        ],
    }

    with open('output_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
