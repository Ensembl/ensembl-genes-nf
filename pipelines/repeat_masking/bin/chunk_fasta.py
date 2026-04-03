#!/usr/bin/env python3
"""
chunk_fasta.py — Split a genome FASTA into fixed-size chunks for
parallelising RepeatMasker, RED, TRF, and DUST.

Each sequence is split into non-overlapping windows of --chunk-size bp.
The last chunk of each sequence may be shorter. Chunk IDs encode the
source sequence and coordinates: seqid_start_end.

Usage:
    chunk_fasta.py --fasta genome.fa --chunk-size 10000000 --outdir chunks/
"""

import argparse
import os
from pathlib import Path


def parse_fasta(path: str):
    """Yield (seqid, sequence) tuples from a FASTA file."""
    seqid, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if seqid is not None:
                    yield seqid, ''.join(parts)
                seqid = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if seqid is not None:
        yield seqid, ''.join(parts)


def write_chunk(seqid: str, start: int, end: int, seq: str, outdir: Path) -> str:
    chunk_id  = f"{seqid}_{start}_{end}"
    out_path  = outdir / f"{chunk_id}.fa"
    with open(out_path, 'w') as fh:
        fh.write(f">{chunk_id}\n")
        chunk_seq = seq[start:end]
        for i in range(0, len(chunk_seq), 60):
            fh.write(chunk_seq[i:i+60] + '\n')
    return str(out_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta',      required=True, help='Input genome FASTA')
    parser.add_argument('--chunk-size', type=int, default=10_000_000,
                        help='Chunk size in bp (default: 10,000,000)')
    parser.add_argument('--outdir',     required=True, help='Output directory for chunk FASTAs')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    total_chunks = 0
    for seqid, seq in parse_fasta(args.fasta):
        seq_len = len(seq)
        for start in range(0, seq_len, args.chunk_size):
            end = min(start + args.chunk_size, seq_len)
            write_chunk(seqid, start, end, seq, outdir)
            total_chunks += 1

    print(f"[chunk_fasta] Wrote {total_chunks} chunks to {outdir}", flush=True)


if __name__ == '__main__':
    main()
