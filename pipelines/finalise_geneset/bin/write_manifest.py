#!/usr/bin/env python3
"""Write output_manifest.json for HiveRunNextflow dataflow."""

import argparse
import json
import os
from datetime import datetime, timezone


def main():
    ap = argparse.ArgumentParser(
        description='Write output_manifest.json for HiveRunNextflow dataflow.'
    )
    ap.add_argument('--pipeline', default='finalise_geneset')
    ap.add_argument('--version',  default='1.0.0')
    ap.add_argument('--outdir',   required=True,
                    help='Pipeline output directory (used to construct the output path)')
    ap.add_argument('--gff3',     required=True,
                    help='Filename (basename) of the GFF3 produced by SELECT_CANONICAL')
    args = ap.parse_args()

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs': [
            {
                'type': 'final_geneset',
                'path': os.path.join(args.outdir, 'finalise_geneset', args.gff3),
            }
        ],
    }

    with open('output_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
