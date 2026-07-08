#!/usr/bin/env python3
"""
Convert RiboTIE transcript-coordinate GTF to genomic coordinates.

RiboTIE outputs predictions in transcript space when run against a
.Aligned.toTranscriptome.out.bam alignment.  This script uses a Gencode
annotation GTF to lift every feature back to genome coordinates.

Multi-exon features are split into one GTF row per genomic block so the
output is compatible with the group_gff() / ribotie_gtf parser in
translon_db_standardise.py.  Output coordinates follow GTF convention
(1-based, inclusive end).

The output filename is the input name with the transcript-coord suffix
stripped (default: '.Aligned.toTranscriptome.out'), so that
collect_ribotie() in build_translon_manifest.py will pick up the genomic
files automatically when they are placed in the deliverables/ tree.

Usage
-----
    python ribotie_tx2genome.py \\
        --gtf    gencode.v47.annotation.gtf.gz \\
        --outdir RiboTIE_results/deliverables/annotated \\
        --input  RiboTIE_results/deliverables/annotated/db_SRR*.Aligned.toTranscriptome.out.annotated.out.gtf

    # Convert both annotated and novel in one call:
    python ribotie_tx2genome.py \\
        --gtf    gencode.v47.annotation.gtf.gz \\
        --outdir RiboTIE_results/deliverables \\
        --input  RiboTIE_results/deliverables/annotated/db_*.Aligned.toTranscriptome.out.annotated.out.gtf \\
                 RiboTIE_results/deliverables/novel/db_*.Aligned.toTranscriptome.out.novel.out.gtf
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def open_maybe_gz(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_gtf_attrs(attr_str: str) -> dict[str, str]:
    """Parse GTF/GFF attribute string, returning key→value dict."""
    result: dict[str, str] = {}
    for m in re.finditer(r'(\w+)\s+"([^"]+)"', attr_str):
        result[m.group(1)] = m.group(2)
    return result


# ---------------------------------------------------------------------------
# Load exon map from Gencode GTF
# ---------------------------------------------------------------------------

def load_exon_map(
    gtf_path: Path,
) -> dict[str, tuple[str, str, list[tuple[int, int]]]]:
    """Parse Gencode GTF and return a transcript → (chrom, strand, exons) map.

    Exons are returned as a list of (start, end) in 0-based half-open
    coordinates, sorted by ascending genomic start.

    Both versioned ('ENST….N') and stripped ('ENST….') IDs are stored so
    lookups succeed regardless of whether the RiboTIE output uses versions.
    """
    # tx_id → {"chrom": str, "strand": str, "exons": list[(int, int)]}
    raw: dict[str, dict] = {}

    print(f"[info] Loading exon map from {gtf_path} ...", file=sys.stderr, flush=True)
    n_exons = 0

    with open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attrs(fields[8])
            tx_id = attrs.get("transcript_id", "")
            if not tx_id:
                continue
            start0 = int(fields[3]) - 1   # GTF 1-based → 0-based
            end0   = int(fields[4])        # GTF inclusive end == half-open end
            if tx_id not in raw:
                raw[tx_id] = {"chrom": fields[0], "strand": fields[6], "exons": []}
            raw[tx_id]["exons"].append((start0, end0))
            n_exons += 1

    # Build final map; also register version-stripped IDs as aliases
    exon_map: dict[str, tuple[str, str, list[tuple[int, int]]]] = {}
    for tx_id, data in raw.items():
        entry = (data["chrom"], data["strand"], sorted(data["exons"]))
        exon_map[tx_id] = entry
        stripped = re.sub(r"\.\d+$", "", tx_id)
        if stripped != tx_id and stripped not in exon_map:
            exon_map[stripped] = entry

    print(
        f"[info] {n_exons:,} exon records → {len(raw):,} transcripts indexed",
        file=sys.stderr, flush=True,
    )
    return exon_map


# ---------------------------------------------------------------------------
# Coordinate conversion
# ---------------------------------------------------------------------------

def tx_interval_to_genomic(
    tx_start0: int,
    tx_end0: int,
    strand: str,
    exons_asc: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Map a 0-based half-open transcript interval to 0-based half-open genomic blocks.

    Parameters
    ----------
    tx_start0, tx_end0 : int
        Feature coordinates in transcript space (0-based half-open).
    strand : str
        Genomic strand of the transcript ('+' or '-').
    exons_asc : list of (int, int)
        Exon intervals (0-based half-open) sorted by ascending genomic start.

    Returns
    -------
    list of (int, int)
        Genomic blocks covering the feature, sorted by ascending genomic start.
        Empty list if the feature does not map onto any exon.

    Notes
    -----
    For '+' strand transcripts, transcript position 0 corresponds to the
    lowest genomic coordinate of the first exon.
    For '-' strand transcripts, transcript position 0 corresponds to the
    highest genomic coordinate of the last exon (in ascending order).
    """
    # Walk exons in transcript order
    tx_ordered = exons_asc if strand == "+" else list(reversed(exons_asc))

    blocks: list[tuple[int, int]] = []
    tx_cursor = 0

    for g_start, g_end in tx_ordered:
        exon_len = g_end - g_start
        exon_tx_start = tx_cursor
        exon_tx_end   = tx_cursor + exon_len

        # Overlap of requested feature with this exon's transcript span
        ov_start = max(tx_start0, exon_tx_start)
        ov_end   = min(tx_end0,   exon_tx_end)

        if ov_start < ov_end:
            rel_s = ov_start - exon_tx_start  # offset into exon from tx start
            rel_e = ov_end   - exon_tx_start

            if strand == "+":
                blocks.append((g_start + rel_s, g_start + rel_e))
            else:
                # For '-' strand, transcript runs from g_end → g_start
                # rel=0 maps to g_end (exclusive), rel=exon_len maps to g_start
                blocks.append((g_end - rel_e, g_end - rel_s))

        tx_cursor = exon_tx_end
        if tx_cursor >= tx_end0:
            break

    return sorted(blocks)


# ---------------------------------------------------------------------------
# Convert a single file
# ---------------------------------------------------------------------------

FEATURE_TYPES = {"CDS", "exon", "translon", "orf"}


def convert_file(
    input_path: Path,
    output_path: Path,
    exon_map: dict[str, tuple[str, str, list[tuple[int, int]]]],
) -> tuple[int, int, int]:
    """Convert one transcript-coord GTF to genomic GTF.

    Returns (n_features_read, n_features_converted, n_skipped).
    n_features_converted counts input features (not output blocks).
    """
    n_read = n_converted = n_skipped = 0

    with open_maybe_gz(input_path) as in_fh, output_path.open("w") as out_fh:
        for line in in_fh:
            if line.startswith("#") or not line.strip():
                out_fh.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                out_fh.write(line)
                continue

            feature = fields[2]
            if feature not in FEATURE_TYPES:
                out_fh.write(line)
                continue

            n_read += 1

            tx_id = fields[0].strip()
            entry = exon_map.get(tx_id) or exon_map.get(re.sub(r"\.\d+$", "", tx_id))
            if entry is None:
                n_skipped += 1
                continue

            chrom, strand, exons_asc = entry

            try:
                tx_start0 = int(fields[3]) - 1
                tx_end0   = int(fields[4])
            except ValueError:
                n_skipped += 1
                continue

            blocks = tx_interval_to_genomic(tx_start0, tx_end0, strand, exons_asc)
            if not blocks:
                n_skipped += 1
                continue

            score = fields[5]
            frame = fields[7]
            attrs = fields[8]

            for g_start0, g_end0 in blocks:
                out_fh.write(
                    f"{chrom}\tRiboTIE\t{feature}\t{g_start0 + 1}\t{g_end0}\t"
                    f"{score}\t{strand}\t{frame}\t{attrs}\n"
                )

            n_converted += 1

    return n_read, n_converted, n_skipped


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input", nargs="+", required=True, type=Path, metavar="GTF",
        help="Transcript-coordinate GTF file(s) to convert (glob-expanded by shell).",
    )
    p.add_argument(
        "--gtf", required=True, type=Path,
        help="Gencode annotation GTF used to build the exon map (may be .gz).",
    )
    p.add_argument(
        "--outdir", required=True, type=Path,
        help=(
            "Output directory.  Converted files are written here with the "
            "transcript-coord suffix removed from the filename.  "
            "Sub-directory structure is NOT preserved; if annotated/ and novel/ "
            "files are mixed, place them in separate --outdir calls."
        ),
    )
    p.add_argument(
        "--tx-suffix", default=".Aligned.toTranscriptome.out",
        metavar="SUFFIX",
        help=(
            "Filename component to strip when building the output name "
            "(default: '.Aligned.toTranscriptome.out').  "
            "Example: 'db_SRR1234.Aligned.toTranscriptome.out.annotated.out.gtf' "
            "→ 'db_SRR1234.annotated.out.gtf'."
        ),
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    exon_map = load_exon_map(args.gtf)

    total_read = total_converted = total_skipped = 0

    for input_path in sorted(args.input):
        out_name = input_path.name.replace(args.tx_suffix, "")
        output_path = args.outdir / out_name

        if output_path == input_path:
            print(
                f"[error] Input and output paths are identical for {input_path}. "
                "Use a different --outdir.",
                file=sys.stderr,
            )
            sys.exit(1)

        n_read, n_converted, n_skipped = convert_file(
            input_path, output_path, exon_map
        )
        pct = f"{100 * n_converted / n_read:.1f}%" if n_read else "n/a"
        print(
            f"  {input_path.name}\n"
            f"    → {output_path}\n"
            f"    features: {n_read} read | {n_converted} converted ({pct}) | {n_skipped} skipped",
            file=sys.stderr,
        )
        total_read      += n_read
        total_converted += n_converted
        total_skipped   += n_skipped

    print(
        f"\n[done] {len(args.input)} file(s) | "
        f"{total_read} features read | "
        f"{total_converted} converted | "
        f"{total_skipped} skipped (transcript not in GTF or no exon overlap)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
