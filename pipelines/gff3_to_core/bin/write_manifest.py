#!/usr/bin/env python3
"""Write output_manifest.json for HiveRunNextflow dataflow."""
import argparse, json, os
from datetime import datetime, timezone

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pipeline', default='gff3_to_core')
    ap.add_argument('--version',  default='1.0.0')
    ap.add_argument('--outdir',   required=True)
    ap.add_argument('--stats-json', required=True)
    ap.add_argument('--dbname',   required=True)
    args = ap.parse_args()

    manifest = {
        'pipeline': args.pipeline,
        'version': args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs': [
            {
                'type': 'core_db',
                'path': args.stats_json,
                'meta': {'dbname': args.dbname},
            }
        ],
    }
    with open('output_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
