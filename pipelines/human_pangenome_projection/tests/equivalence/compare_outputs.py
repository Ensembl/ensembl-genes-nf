#!/usr/bin/env python3
"""Compare two pipeline output directories for byte-level equivalence.

Used to prove that the modular / Nextflow-orchestrated pipeline produces output
identical to the in-memory monolith. Two artifact families are compared:

  * The mapped GFF3 — required to be byte-identical (its only header line is
    ``##gff-version 3``; the body is coordinate-sorted by the writer).
  * The stats / audit JSON and the removed-genes TSV — identical after stripping
    run-specific volatile fields (wall-clock ``timestamp`` and ``*_time_seconds``
    timings), which carry no biological meaning.

Exit code 0 means equivalent; non-zero prints the first differences found.

Usage:
    compare_outputs.py REF_DIR TEST_DIR \
        --gff chm13.chr20.mapped.gff3 \
        --stats chr20.stats.json
"""

import argparse
import json
import sys
from pathlib import Path

# Keys removed before comparing any JSON (recursively): wall-clock + timings.
_VOLATILE_KEYS = {
    "timestamp",
    "total_time_seconds",
    "synteny_time_seconds",
    "mapping_time_seconds",
    "validation_time_seconds",
}


def _strip_volatile(obj):
    if isinstance(obj, dict):
        return {k: _strip_volatile(v) for k, v in obj.items() if k not in _VOLATILE_KEYS}
    if isinstance(obj, list):
        return [_strip_volatile(x) for x in obj]
    return obj


def _read_text(path: Path) -> str:
    return path.read_text() if path.exists() else None


def compare_gff(ref: Path, test: Path) -> list:
    diffs = []
    a, b = _read_text(ref), _read_text(test)
    if a is None or b is None:
        return [f"missing GFF: ref_exists={a is not None} test_exists={b is not None}"]
    if a == b:
        return []
    a_lines, b_lines = a.splitlines(), b.splitlines()
    if len(a_lines) != len(b_lines):
        diffs.append(f"GFF line count differs: ref={len(a_lines)} test={len(b_lines)}")
    for i, (la, lb) in enumerate(zip(a_lines, b_lines), 1):
        if la != lb:
            diffs.append(f"GFF first diff at line {i}:\n  ref:  {la}\n  test: {lb}")
            break
    return diffs


def compare_json(ref: Path, test: Path) -> list:
    if not ref.exists() and not test.exists():
        return []
    if not ref.exists() or not test.exists():
        return [f"missing JSON: ref_exists={ref.exists()} test_exists={test.exists()} ({ref.name})"]
    a = _strip_volatile(json.loads(ref.read_text()))
    b = _strip_volatile(json.loads(test.read_text()))
    if a == b:
        return []
    # Locate first differing top-level key for a useful message.
    diffs = [f"JSON differs after stripping volatile keys: {ref.name}"]
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if a.get(key) != b.get(key):
                diffs.append(f"  key '{key}': ref={a.get(key)!r:.200} test={b.get(key)!r:.200}")
    return diffs


def compare_text(ref: Path, test: Path) -> list:
    if not ref.exists() and not test.exists():
        return []
    a, b = _read_text(ref), _read_text(test)
    if a == b:
        return []
    return [f"text differs: {ref.name}"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ref_dir", type=Path)
    ap.add_argument("test_dir", type=Path)
    ap.add_argument("--gff", default="chm13.chr20.mapped.gff3")
    ap.add_argument("--stats", default="chr20.stats.json")
    args = ap.parse_args()

    stats_stem = args.stats[: -len(".json")] if args.stats.endswith(".json") else args.stats

    checks = []
    checks += [("GFF3", compare_gff(args.ref_dir / args.gff, args.test_dir / args.gff))]
    checks += [("stats.json", compare_json(args.ref_dir / args.stats, args.test_dir / args.stats))]
    checks += [("audit.json", compare_json(
        args.ref_dir / f"{stats_stem}.audit.json", args.test_dir / f"{stats_stem}.audit.json"))]
    checks += [("removed.tsv", compare_text(
        args.ref_dir / f"{stats_stem}.removed.tsv", args.test_dir / f"{stats_stem}.removed.tsv"))]

    failed = False
    for name, diffs in checks:
        if diffs:
            failed = True
            print(f"[DIFF] {name}")
            for d in diffs:
                print(f"       {d}")
        else:
            print(f"[ OK ] {name}")

    if failed:
        print("\nRESULT: outputs DIFFER")
        return 1
    print("\nRESULT: outputs are EQUIVALENT (byte-identical modulo timings)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
