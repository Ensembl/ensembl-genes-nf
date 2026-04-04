#!/usr/bin/env python3
"""write_manifest.py for repeat_masking pipeline."""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pipeline',          default='repeat_masking')
    parser.add_argument('--version',           default='1.0.0')
    parser.add_argument('--outdir',            required=True)
    parser.add_argument('--softmasked-fasta',  required=True)
    parser.add_argument('--repeat-gff3',       required=True)
    args = parser.parse_args()

    # Build outdir-relative paths so the manifest is stable after the
    # Nextflow work directory is cleaned.  BEDTOOLS_MASKFASTA publishes
    # its *.softmasked.fa to ${outdir}/genome/ via publishDir.
    outputs = []
    for path_str, out_type, subdir in [
        (args.softmasked_fasta, 'softmasked_fasta', 'genome'),
        (args.repeat_gff3,      'repeat_gff3',      'repeats'),
    ]:
        p = Path(path_str)
        published = Path(args.outdir) / subdir / p.name
        outputs.append({
            'type': out_type,
            'path': str(published),
            'meta': {'id': p.stem.split('.')[0]}
        })

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'outputs':      outputs,
    }

    with open('output_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)

    print(f"[write_manifest] Written {len(outputs)} outputs", flush=True)


if __name__ == '__main__':
    main()
