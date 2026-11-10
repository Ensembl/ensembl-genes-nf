#!/usr/bin/env python3
"""Record legacy features and direct alignment support in one streamed scan."""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import shutil
import sqlite3
import subprocess
from contextlib import contextmanager
from pathlib import Path

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
LEGACY_EXCLUDED = 4 | 256 | 512 | 1024 | 2048
SUPPORT_EXCLUDED = 4 | 256 | 2048
SCHEMA = "direct-alignment-support-v1"


@contextmanager
def gzip_text(path: Path):
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as packed:
            with io.TextIOWrapper(packed, encoding="utf-8", newline="") as handle:
                yield handle


@contextmanager
def sam_lines(bam: Path | None, sam: Path | None):
    if sam is not None:
        opener = gzip.open if sam.suffix == ".gz" else Path.open
        with opener(sam, "rt", encoding="utf-8") as handle:
            yield handle
        return
    if shutil.which("samtools") is None:
        raise SystemExit("samtools is required for BAM extraction")
    proc = subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE,
                            text=True, encoding="utf-8", bufsize=1024 * 1024)
    assert proc.stdout is not None
    try:
        yield proc.stdout
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            raise SystemExit(f"samtools view failed: {bam}")


def geometry(left: int, cigar: str):
    parts = CIGAR.findall(cigar)
    if not parts or "".join(f"{n}{op}" for n, op in parts) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar}")
    right, introns, aligned = left, [], 0
    for number, op in parts:
        size = int(number)
        if op == "N":
            introns.append((right, right + size))
            right += size
        elif op in "MD=X":
            right += size
        if op in "M=X":
            aligned += size
    chain = ";".join(f"{start}-{end}" for start, end in introns)
    left_clip = int(parts[0][0]) if parts[0][1] in "SH" else 0
    right_clip = int(parts[-1][0]) if parts[-1][1] in "SH" else 0
    return right, introns, chain, aligned, left_clip, right_clip


def table(path: Path, fields, rows):
    writer_context = gzip_text(path) if path.suffix == ".gz" else path.open("w", newline="")
    with writer_context as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", type=Path)
    source.add_argument("--sam", type=Path)
    for name in ("read_features", "five_prime", "three_prime", "junctions", "summary", "support_dir", "run_id"):
        parser.add_argument("--" + name.replace("_", "-"), required=True, type=Path if name != "run_id" else str)
    for name in ("assembly_id", "sample_id", "library_id", "read_set_id", "technology"):
        parser.add_argument("--" + name.replace("_", "-"), default="unknown")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--skip-support", action="store_true")
    args = parser.parse_args()
    if args.min_mapq < 0:
        parser.error("--min-mapq must be nonnegative")
    args.support_dir.mkdir(parents=True, exist_ok=True)
    db_path = args.support_dir / "observations.sqlite"
    db = sqlite3.connect(db_path)
    db.execute("CREATE TABLE observations (contig TEXT, chain TEXT, orientation TEXT, left_end INTEGER, right_end INTEGER, read_id TEXT, mapq INTEGER, aligned_bases INTEGER, left_clip INTEGER, right_clip INTEGER, nm INTEGER)")
    db.execute("CREATE TABLE junctions (contig TEXT, start INTEGER, end INTEGER, orientation TEXT, chain TEXT, position INTEGER, read_id TEXT)")
    db.execute("CREATE TABLE legacy_junctions (contig TEXT, start INTEGER, end INTEGER, orientation TEXT)")
    keys = ("records_scanned", "primary_mapped", "mapq_passing", "legacy_retained", "chain_records", "junction_records", "end_records", "duplicate_records", "qc_fail_records", "secondary_records", "supplementary_records", "unmapped_records")
    flow = dict.fromkeys(keys, 0)
    with sam_lines(args.bam, args.sam) as records, gzip_text(args.read_features) as features, args.five_prime.open("w") as five, args.three_prime.open("w") as three:
        features.write("read_id\tcontig\tstrand\tstart\tend\taligned_bases\tjunction_count\tcigar\n")
        for line in records:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise SystemExit("malformed SAM alignment record")
            flow["records_scanned"] += 1
            read_id, flag, contig = fields[0], int(fields[1]), fields[2]
            mapq, cigar = int(fields[4]), fields[5]
            for bit, key in ((4, "unmapped_records"), (256, "secondary_records"), (512, "qc_fail_records"), (1024, "duplicate_records"), (2048, "supplementary_records")):
                if flag & bit:
                    flow[key] += 1
            if flag & 4:
                continue
            if not flag & (256 | 2048):
                flow["primary_mapped"] += 1
            left = int(fields[3]) - 1
            right, introns, chain, aligned, left_clip, right_clip = geometry(left, cigar)
            orientation = "-" if flag & 16 else "+"
            if not flag & LEGACY_EXCLUDED:
                flow["legacy_retained"] += 1
                five_end = left if orientation == "+" else right - 1
                three_end = right - 1 if orientation == "+" else left
                features.write(f"{read_id}\t{contig}\t{orientation}\t{left}\t{right}\t{aligned}\t{len(introns)}\t{cigar}\n")
                five.write(f"{contig}\t{five_end}\t{five_end + 1}\t{read_id}\t0\t{orientation}\n")
                three.write(f"{contig}\t{three_end}\t{three_end + 1}\t{read_id}\t0\t{orientation}\n")
                db.executemany("INSERT INTO legacy_junctions VALUES (?,?,?,?)", ((contig, start, end, orientation) for start, end in introns))
            if args.skip_support or flag & SUPPORT_EXCLUDED or mapq < args.min_mapq:
                continue
            flow["mapq_passing"] += 1
            flow["end_records"] += 1
            if introns:
                flow["chain_records"] += 1
                flow["junction_records"] += 1
            nm = next((int(tag[5:]) for tag in fields[11:] if tag.startswith("NM:i:")), None)
            db.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?,?,?,?)", (contig, chain, orientation, left, right - 1, read_id, mapq, aligned, left_clip, right_clip, nm))
            db.executemany("INSERT INTO junctions VALUES (?,?,?,?,?,?,?)", ((contig, start, end, orientation, chain, position, read_id) for position, (start, end) in enumerate(introns, 1)))
    db.commit()
    with args.junctions.open("w") as handle:
        for contig, start, end, orientation, count in db.execute("SELECT contig,start,end,orientation,COUNT(*) FROM legacy_junctions GROUP BY contig,start,end,orientation ORDER BY contig,start,end,orientation"):
            handle.write(f"{contig}\t{start}\t{end}\t{contig}:{start}-{end}\t{count}\t{orientation}\n")
    distinct = db.execute("SELECT COUNT(*) FROM (SELECT 1 FROM legacy_junctions GROUP BY contig,start,end,orientation)").fetchone()[0]
    args.summary.write_text(f"run_id\treads\tjunctions\n{args.run_id}\t{flow['legacy_retained']}\t{distinct}\n")
    if args.skip_support:
        db.close()
        db_path.unlink()
        return
    prefix = (SCHEMA, args.run_id, args.assembly_id, args.sample_id, args.library_id, args.read_set_id, args.technology)
    common = ("schema_version", "run_id", "assembly_id", "sample_id", "library_id", "read_set_id", "technology")
    chains = ((*prefix, contig, chain, orientation, chain.count(";") + 1 if chain else 0, count, reads, min_left, max_right, min_mapq, max_mapq)
              for contig, chain, orientation, count, reads, min_left, max_right, min_mapq, max_mapq in db.execute("SELECT contig,chain,orientation,COUNT(*),COUNT(DISTINCT read_id),MIN(left_end),MAX(right_end),MIN(mapq),MAX(mapq) FROM observations GROUP BY contig,chain,orientation ORDER BY contig,chain,orientation"))
    table(args.support_dir / "chain_observations.tsv.gz", common + ("contig", "intron_chain", "alignment_orientation", "intron_count", "alignment_records", "distinct_read_ids", "min_left_end", "max_right_end", "min_mapq", "max_mapq"), chains)
    junctions = ((*prefix, contig, start, end, orientation, chain, position, count, reads, "not_available")
                 for contig, start, end, orientation, chain, position, count, reads in db.execute("SELECT contig,start,end,orientation,chain,position,COUNT(*),COUNT(DISTINCT read_id) FROM junctions GROUP BY contig,start,end,orientation,chain,position ORDER BY contig,start,end,orientation,chain,position"))
    table(args.support_dir / "junction_observations.tsv.gz", common + ("contig", "intron_start", "intron_end", "alignment_orientation", "intron_chain", "chain_position", "alignment_records", "distinct_read_ids", "splice_motif_status"), junctions)
    ends = ((*prefix, contig, chain, orientation, left, right, count, reads)
            for contig, chain, orientation, left, right, count, reads in db.execute("SELECT contig,chain,orientation,left_end,right_end,COUNT(*),COUNT(DISTINCT read_id) FROM observations GROUP BY contig,chain,orientation,left_end,right_end ORDER BY contig,chain,orientation,left_end,right_end"))
    table(args.support_dir / "end_observations.tsv.gz", common + ("contig", "intron_chain", "alignment_orientation", "left_aligned_end", "right_aligned_end", "alignment_records", "distinct_read_ids"), ends)
    table(args.support_dir / "alignment_filter_flow.tsv", common + ("stage", "alignment_records"), ((*prefix, key, value) for key, value in flow.items()))
    noise = db.execute("SELECT COUNT(*),SUM(left_clip>0),SUM(right_clip>0),SUM(nm IS NOT NULL),MIN(mapq),MAX(mapq),AVG(mapq),AVG(nm) FROM observations").fetchone()
    table(args.support_dir / "alignment_noise.tsv", common + ("retained_records", "left_clipped_records", "right_clipped_records", "nm_available_records", "min_mapq", "max_mapq", "mean_mapq", "mean_nm"), [(*prefix, *(value if value is not None else "NA" for value in noise))])
    metadata = {"schema": SCHEMA, "run_id": args.run_id, "assembly_id": args.assembly_id, "sample_id": args.sample_id, "library_id": args.library_id, "read_set_id": args.read_set_id, "technology": args.technology, "coordinate_convention": "zero_based_half_open_introns;aligned_end_base_zero_based", "support_filter": f"-F {SUPPORT_EXCLUDED} -q {args.min_mapq}", "legacy_filter": f"-F {LEGACY_EXCLUDED}", "full_length_status": "unknown", "molecule_status": "not_available", "alignment_orientation_is_transcript_strand": False}
    (args.support_dir / "evidence_metadata.json").write_text(json.dumps(metadata, sort_keys=True, indent=2) + "\n")
    db.close()
    db_path.unlink()


if __name__ == "__main__":
    main()
