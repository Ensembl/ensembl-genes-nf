#!/usr/bin/env python3
"""Write output_manifest.json for HiveRunNextflow dataflow."""

import argparse
import json
import os
from datetime import datetime, timezone


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pipeline',    default='load_assembly')
    ap.add_argument('--version',     default='1.0.0')
    ap.add_argument('--outdir',      required=True)
    ap.add_argument('--genome_fasta', required=True)
    ap.add_argument('--synonyms_tsv', required=True)
    ap.add_argument('--meta_json',    required=True)
    args = ap.parse_args()

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs': [
            {
                'type': 'genome_fasta',
                'path': os.path.join(args.outdir, 'genome', args.genome_fasta),
            },
            {
                'type': 'seq_region_synonyms',
                'path': os.path.join(args.outdir, 'genome', args.synonyms_tsv),
            },
            {
                'type': 'assembly_metadata',
                'path': os.path.join(args.outdir, 'genome', args.meta_json),
            },
        ],
    }

    with open('output_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
