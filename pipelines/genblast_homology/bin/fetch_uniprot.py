#!/usr/bin/env python3
"""
Download reviewed UniProt proteins for a taxon using the UniProt REST API.

UniProt REST streaming endpoint:
  https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(taxonomy_id:TAXON)+AND+(reviewed:true)

Usage:
  fetch_uniprot.py --taxon-id 9989 --out rodentia_reviewed.fa
  fetch_uniprot.py --taxon-id 40674 --out mammalia_reviewed.fa  # larger set
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
MAX_RETRIES = 5
RETRY_DELAY = 30  # seconds — UniProt can be slow


def fetch(taxon_id: str, out_path: Path) -> int:
    """Stream UniProt FASTA for a taxon. Returns number of sequences."""
    query  = f"(taxonomy_id:{taxon_id}) AND (reviewed:true)"
    params = {"format": "fasta", "query": query}

    for attempt in range(MAX_RETRIES):
        try:
            print(f"Querying UniProt for taxon {taxon_id} (attempt {attempt+1})...", file=sys.stderr)
            resp = requests.get(UNIPROT_STREAM_URL, params=params, stream=True, timeout=600)
            resp.raise_for_status()

            count = 0
            with open(out_path, "w") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20, decode_unicode=True):
                    fh.write(chunk)
                    count += chunk.count(">")

            print(f"Downloaded ~{count} sequences to {out_path}", file=sys.stderr)
            return count

        except Exception as exc:
            if attempt == MAX_RETRIES - 1:
                raise RuntimeError(
                    f"UniProt download failed after {MAX_RETRIES} attempts: {exc}"
                ) from exc
            print(f"  Retry in {RETRY_DELAY}s: {exc}", file=sys.stderr)
            time.sleep(RETRY_DELAY)

    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxon-id", required=True,  help="NCBI taxonomy ID (integer)")
    parser.add_argument("--out",      required=True,  help="Output FASTA path")
    args = parser.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    n = fetch(args.taxon_id, out)
    if n == 0:
        print(f"WARNING: no sequences downloaded for taxon {args.taxon_id}. "
              f"Check the taxon ID and UniProt availability.", file=sys.stderr)
        # Write an empty file so Nextflow doesn't fail on missing output
        out.touch()


if __name__ == "__main__":
    main()
