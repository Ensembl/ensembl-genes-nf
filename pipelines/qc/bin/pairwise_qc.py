#!/usr/bin/env python3
"""Format-neutral helpers for pairwise transcript-model QC.

The command line is intentionally small: Nextflow modules invoke one
subcommand per task, while parsing, locus assignment, rendering and metric
contracts live here. The canonical table is the boundary between source
formats and every structural comparison.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
import csv
import gzip
import hashlib
import heapq
from itertools import islice
import shutil
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse
from urllib.parse import unquote
from urllib.request import urlopen


MODEL_HEADER = [
    "pair_id", "source", "assembly", "contig", "strand", "transcript_id",
    "original_transcript_id", "original_gene_id", "transcript_start",
    "transcript_end", "exon_count", "exon_coordinates", "intron_chain",
    "five_prime_end", "three_prime_end", "source_format",
    "gene_assignment_type", "locus_id", "locus_assignment_method",
    "model_status", "utr_coordinates", "cds_coordinates",
]

ASSIGNMENT_TYPES = {
    "biological_source_gene", "reference_gene", "synthetic_locus", "unassigned"
}
VALID_FORMATS = {"gff3", "gtf", "bed12", "auto"}


def open_text(path: str | Path):
    path = str(path)
    with open(path, "rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8") if compressed or path.endswith(".gz") else open(path, encoding="utf-8")


def attrs(text: str) -> dict[str, str]:
    """Parse the common GTF and GFF3 attribute forms."""
    result: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            match = re.match(r"([^\s]+)\s+(?:\"([^\"]*)\"|([^\s]+))$", item)
            if not match:
                continue
            key, value = match.group(1), match.group(2) or match.group(3)
        result[key.strip()] = unquote(value.strip().strip('"'))
    return result


def fmt_coord(start: int, end: int) -> str:
    return f"{start}-{end}"


def parse_coord(value: str) -> tuple[int, int] | None:
    try:
        start, end = value.split("-", 1)
        return int(start), int(end)
    except (AttributeError, ValueError):
        return None


def intervals(row: dict[str, str], field: str = "exon_coordinates") -> list[tuple[int, int]]:
    if field == "exon_coordinates" and "_parsed_exons" in row:
        return row["_parsed_exons"]
    parsed = sorted(coord for coord in (parse_coord(x) for x in row.get(field, "").split(";")) if coord)
    if field == "exon_coordinates":
        row["_parsed_exons"] = parsed
    return parsed


def chain(exons: list[tuple[int, int]]) -> str:
    ordered = sorted(exons)
    return ";".join(f"{right}-{next_left}" for (_, right), (next_left, _) in zip(ordered, ordered[1:]))


def junctions(row: dict[str, str]) -> set[str]:
    if "_parsed_junctions" in row:
        return row["_parsed_junctions"]
    parsed = {x for x in row.get("intron_chain", "").split(";") if x}
    row["_parsed_junctions"] = parsed
    return parsed


def infer_format(path: str, requested: str | None = None) -> str:
    if requested and requested.lower() not in {"auto", ""}:
        return requested.lower()
    name = str(path).lower()
    if name.endswith((".bed", ".bed12", ".bed.gz", ".bed12.gz")):
        return "bed12"
    if name.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    return "gff3"


def valid_gene_id(value: str | None) -> bool:
    return bool(value and value.strip() and value.strip() not in {".", "-", "NA", "null"}
                and not value.startswith("synthetic_gene:"))


def empty_assignment(row: dict[str, str]) -> dict[str, str]:
    row = dict(row)
    row.setdefault("transcript_id", row.get("original_transcript_id", ""))
    row.setdefault("original_transcript_id", row.get("transcript_id", ""))
    row.setdefault("gene_assignment_type", "unassigned")
    row.setdefault("locus_id", "")
    row.setdefault("locus_assignment_method", "")
    return row


def safe_int(value: str | None, default: int = 10**18) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def coordinate_span(row: dict[str, str]) -> tuple[int, int] | None:
    """Return a model's genomic span when both coordinates are numeric."""
    try:
        start, end = int(row.get("transcript_start", "")), int(row.get("transcript_end", ""))
    except (TypeError, ValueError):
        return None
    return (start, end) if start <= end else None


def build_interval_index(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], tuple[list[int], list[tuple[int, int, int]]]]:
    """Index transcript spans by assembly, contig and strand for overlap queries."""
    buckets: defaultdict[tuple[str, str, str], list[tuple[int, int, int]]] = defaultdict(list)
    for index, row in enumerate(rows):
        span = coordinate_span(row)
        if span is None:
            continue
        key = (row.get("assembly", ""), row.get("contig", ""), row.get("strand", "."))
        buckets[key].append((span[0], span[1], index))

    indexed = {}
    for key, entries in buckets.items():
        entries.sort(key=lambda item: (item[0], item[1], item[2]))
        indexed[key] = ([item[0] for item in entries], entries)
    return indexed


def indexed_candidates(row: dict[str, str], index: dict[tuple[str, str, str], tuple[list[int], list[tuple[int, int, int]]]]):
    """Yield row indexes whose transcript spans overlap the query span."""
    span = coordinate_span(row)
    if span is None:
        return
    key = (row.get("assembly", ""), row.get("contig", ""), row.get("strand", "."))
    bucket = index.get(key)
    if bucket is None:
        return
    starts, entries = bucket
    # Do not materialise a slice for every query. Large BED12 inputs can
    # contain millions of models, and the slice itself becomes a substantial
    # allocation even when only a few entries overlap the query span.
    for start, end, row_index in islice(entries, bisect_right(starts, span[1])):
        if end >= span[0]:
            yield row_index


def parse_models(path: str, source: str, pair_id: str, assembly: str,
                 requested: str | None = None) -> tuple[list[dict[str, str]], list[str]]:
    fmt = infer_format(path, requested)
    if fmt not in VALID_FORMATS - {"auto"}:
        raise ValueError(f"unsupported source format: {fmt}")

    records: dict[str, dict] = {}
    warnings: list[str] = []
    if fmt == "bed12":
        with open_text(path) as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    warnings.append(f"bed12_line_{line_no}: fewer than 12 columns")
                    continue
                try:
                    chrom, start0, end, name, _score, strand = fields[:6]
                    start0, end, block_count = int(start0), int(end), int(fields[9])
                    sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                    offsets = [int(x) for x in fields[11].rstrip(",").split(",") if x]
                    if len(sizes) != block_count or len(offsets) != block_count:
                        raise ValueError("blockCount does not match block lists")
                    exons = [(start0 + offset + 1, start0 + offset + size)
                             for offset, size in zip(offsets, sizes)]
                    if not exons or min(exons)[0] < 1 or max(exons)[1] > end:
                        raise ValueError("invalid BED12 block coordinates")
                    if any(size <= 0 or offset < 0 for size, offset in zip(sizes, offsets)):
                        raise ValueError("BED12 blocks must have positive sizes and non-negative offsets")
                except (ValueError, IndexError) as exc:
                    warnings.append(f"bed12_line_{line_no}: {exc}")
                    continue
                tx = name or f"bed12_tx_{line_no}"
                records[tx] = {
                    "gene": "", "tx": tx, "contig": chrom, "strand": strand or ".",
                    "exons": sorted(exons), "utr": [], "cds": [], "status": "valid",
                }
    else:
        tx_types = {
            "transcript", "mrna", "lnc_rna", "ncrna", "mirna", "snrna", "snorna",
            "pseudogenic_transcript", "antisense_rna", "trna", "rrna", "scrna",
            "guide_rna", "vault_rna", "y_rna",
        }
        pending_exons: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
        pending_utr: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
        pending_cds: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
        tx_gene: dict[str, str] = {}
        tx_meta: dict[str, dict] = {}
        gene_ids: dict[str, str] = {}

        with open_text(path) as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    warnings.append(f"annotation_line_{line_no}: expected 9 columns")
                    continue
                contig, _src, feature, start, end, _score, strand, _phase, raw_attrs = fields
                try:
                    start_i, end_i = int(start), int(end)
                    if start_i < 1 or end_i < start_i:
                        raise ValueError("invalid coordinates")
                except ValueError as exc:
                    warnings.append(f"annotation_line_{line_no}: {exc}")
                    continue
                a = attrs(raw_attrs)
                feature_l = feature.lower()
                feature_id = a.get("ID") or a.get("gene_id") or a.get("gene")
                if feature_l == "gene" and feature_id:
                    gene_ids[feature_id] = a.get("gene_id") or feature_id
                    continue
                if feature_l in tx_types:
                    tx = a.get("transcript_id") or a.get("ID") or a.get("transcript")
                    parent = (a.get("Parent") or "").split(",")[0]
                    gene = (a.get("gene_id") or a.get("gene") or a.get("gene_name")
                            or gene_ids.get(parent) or (parent if fmt == "gff3" else ""))
                    if tx:
                        tx_gene[tx] = gene
                        tx_meta[tx] = {
                            "contig": contig, "strand": strand or ".", "start": start_i,
                            "end": end_i, "status": "valid",
                        }
                elif feature_l == "exon":
                    parents = (a.get("Parent") or a.get("transcript_id") or a.get("transcript") or "").split(",")
                    for tx in filter(None, parents):
                        pending_exons[tx].append((start_i, end_i))
                        if tx not in tx_meta:
                            tx_meta[tx] = {"contig": contig, "strand": strand or "."}
                        tx_gene.setdefault(tx, a.get("gene_id") or a.get("gene") or "")
                elif feature_l in {"utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr"}:
                    for tx in filter(None, (a.get("Parent") or a.get("transcript_id") or "").split(",")):
                        pending_utr[tx].append((start_i, end_i))
                elif feature_l == "cds":
                    for tx in filter(None, (a.get("Parent") or a.get("transcript_id") or "").split(",")):
                        pending_cds[tx].append((start_i, end_i))

        for tx, gene in tx_gene.items():
            exons = sorted(pending_exons.get(tx, []))
            meta = tx_meta.get(tx, {})
            if not exons and meta.get("start") is not None:
                exons = [(meta["start"], meta["end"])]
                status = "no_exon_records"
            elif not exons:
                status = "missing_exons"
            else:
                status = meta.get("status", "valid")
            start_i = min(x[0] for x in exons) if exons else meta.get("start", "")
            end_i = max(x[1] for x in exons) if exons else meta.get("end", "")
            records[tx] = {
                "gene": gene if valid_gene_id(gene) else "", "tx": tx,
                "contig": meta.get("contig", ""), "strand": meta.get("strand", "."),
                "exons": exons, "utr": sorted(pending_utr.get(tx, [])),
                "cds": sorted(pending_cds.get(tx, [])), "status": status,
                "start": start_i, "end": end_i,
            }

    rows: list[dict[str, str]] = []
    for rec in records.values():
        exons = sorted(rec.get("exons", []))
        start_i = rec.get("start") or (min(x[0] for x in exons) if exons else "")
        end_i = rec.get("end") or (max(x[1] for x in exons) if exons else "")
        five, three = (end_i, start_i) if rec.get("strand") == "-" else (start_i, end_i)
        gene = rec.get("gene", "")
        rows.append(empty_assignment({
            "pair_id": pair_id, "source": source, "assembly": assembly,
            "contig": rec.get("contig", ""), "strand": rec.get("strand", "."),
            "transcript_id": rec["tx"], "original_transcript_id": rec["tx"],
            "original_gene_id": gene, "transcript_start": str(start_i),
            "transcript_end": str(end_i), "exon_count": str(len(exons)),
            "exon_coordinates": ";".join(fmt_coord(*x) for x in exons),
            "intron_chain": chain(exons), "five_prime_end": str(five),
            "three_prime_end": str(three), "source_format": fmt,
            "gene_assignment_type": "biological_source_gene" if valid_gene_id(gene) else "unassigned",
            "locus_id": "", "locus_assignment_method": "",
            "model_status": rec.get("status", "valid"),
            "utr_coordinates": ";".join(fmt_coord(*x) for x in rec.get("utr", [])),
            "cds_coordinates": ";".join(fmt_coord(*x) for x in rec.get("cds", [])),
        }))
    return sorted(rows, key=lambda row: (row["contig"], safe_int(row["transcript_start"]), row["transcript_id"])), warnings


def iter_tsv(path: str | Path):
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return
    with open(path, newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    return list(iter_tsv(path))


def write_tsv(path: str | Path, header: list[str], rows: list[dict | list]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([row.get(column, "") for column in header] if isinstance(row, dict) else row)


def stage_input(value: str | None, destination: str | Path, cache_dir: str | Path) -> None:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    if not value or value.lower() in {"na", "null"}:
        destination.touch()
        return
    if re.match(r"(?i)^(https?|ftp)://", value):
        cache = Path(cache_dir)
        cache.mkdir(parents=True, exist_ok=True)
        suffix = Path(urlparse(value).path).suffix
        cached = cache / f"{hashlib.sha256(value.encode('utf-8')).hexdigest()}{suffix}"
        if not cached.is_file() or cached.stat().st_size == 0:
            with urlopen(value) as response, open(cached, "wb") as handle:
                shutil.copyfileobj(response, handle)
        shutil.copyfile(cached, destination)
        return
    source = Path(value)
    if not source.is_file():
        raise SystemExit(f"input does not exist: {value}")
    shutil.copyfile(source, destination)


def ordered_exons(row: dict[str, str]) -> list[tuple[int, int]]:
    exons = intervals(row)
    return list(reversed(exons)) if row.get("strand") == "-" else exons


def feature_signature(row: dict[str, str]) -> dict[str, object]:
    exons = intervals(row)
    ordered = ordered_exons(row)
    contig = row.get("contig", "")
    strand = row.get("strand", ".")
    exon_keys = [f"{contig}|{strand}|{start}-{end}" for start, end in exons]
    genomic_junctions = [
        f"{right}-{next_left}"
        for (_, right), (next_left, _) in zip(exons, exons[1:])
    ]
    chain_keys = list(genomic_junctions if strand != "-" else reversed(genomic_junctions))
    chain_key = f"{contig}|{strand}|{';'.join(chain_keys)}" if chain_keys else ""
    model_key = f"{contig}|{strand}|{';'.join(f'{start}-{end}' for start, end in exons)}"
    five = ordered[0] if ordered else None
    three = ordered[-1] if ordered else None
    five_key = f"5p|{contig}|{strand}|{five[0]}-{five[1]}" if five else ""
    three_key = f"3p|{contig}|{strand}|{three[0]}-{three[1]}" if three else ""
    return {
        "model": model_key if exons else "",
        "exons": exon_keys,
        "introns": [f"{contig}|{strand}|{junction}" for junction in genomic_junctions],
        "intron_chain": chain_key,
        "five_prime_exon": five_key,
        "three_prime_exon": three_key,
    }


def add_example(mapping: defaultdict[str, list[str]], key: str, value: str, limit: int = 5) -> None:
    if key and len(mapping[key]) < limit and value not in mapping[key]:
        mapping[key].append(value)


def build_feature_index(path: str | Path) -> dict[str, object]:
    index = {
        "features": defaultdict(Counter),
        "examples": defaultdict(list),
        "chain_terminals": defaultdict(set),
        "model_count": 0,
    }
    for row in iter_tsv(path):
        signature = feature_signature(row)
        transcript_id = row.get("transcript_id", row.get("original_transcript_id", ""))
        index["model_count"] += 1
        for feature_type, key in (("model", signature["model"]), ("intron_chain", signature["intron_chain"])):
            if key:
                index["features"][feature_type][key] += 1
                add_example(index["examples"], key, transcript_id)
        for feature_type in ("exons", "introns"):
            for key in signature[feature_type]:
                index["features"][feature_type[:-1]][key] += 1
                add_example(index["examples"], key, transcript_id)
        for feature_type, key in (("terminal_exon", signature["five_prime_exon"]), ("terminal_exon", signature["three_prime_exon"])):
            if key:
                index["features"][feature_type][key] += 1
                add_example(index["examples"], key, transcript_id)
        if signature["intron_chain"]:
            index["chain_terminals"][signature["intron_chain"]].add((signature["five_prime_exon"], signature["three_prime_exon"]))
    return index


def classify_model_difference(signature: dict[str, object], other: dict[str, object]) -> tuple[str, int, int, int, int]:
    model_key = signature["model"]
    chain_key = signature["intron_chain"]
    exact_models = other["features"]["model"].get(model_key, 0) if model_key else 0
    same_chain_models = other["features"]["intron_chain"].get(chain_key, 0) if chain_key else 0
    shared_introns = sum(1 for key in signature["introns"] if other["features"]["intron"].get(key, 0))
    shared_exons = sum(1 for key in signature["exons"] if other["features"]["exon"].get(key, 0))
    if exact_models:
        classification = "exact_model"
    elif same_chain_models:
        terminals = other["chain_terminals"].get(chain_key, set())
        same_five = any(pair[0] == signature["five_prime_exon"] for pair in terminals)
        same_three = any(pair[1] == signature["three_prime_exon"] for pair in terminals)
        if same_five and not same_three:
            classification = "same_intron_chain_3prime_terminal_difference"
        elif same_three and not same_five:
            classification = "same_intron_chain_5prime_terminal_difference"
        else:
            classification = "same_intron_chain_terminal_difference"
    elif shared_introns:
        classification = "shared_introns_different_chain"
    elif shared_exons:
        classification = "shared_exons_only"
    else:
        classification = "no_shared_structure"
    return classification, exact_models, same_chain_models, shared_introns, shared_exons


def model_presence_rows(path: str | Path, other: dict[str, object], pair_id: str, source: str,
                        evidence: dict[str, dict[str, str]]):
    for row in iter_tsv(path):
        signature = feature_signature(row)
        tx = row.get("transcript_id", row.get("original_transcript_id", ""))
        ev = evidence.get(tx, {})
        classification, exact, same_chain, shared_introns, shared_exons = classify_model_difference(signature, other)
        yield [
            pair_id, source, tx, signature["model"], ev.get("full_length_support", ""),
            "true" if ev else "false", str(exact), str(same_chain), str(shared_introns),
            str(shared_exons), classification,
        ]


def write_model_presence(path: str | Path, pair_id: str, source_a: str | Path, source_b: str | Path,
                         index_a: dict[str, object], index_b: dict[str, object], evidence_a: dict[str, dict[str, str]],
                         evidence_b: dict[str, dict[str, str]]) -> None:
    header = [
        "pair_id", "source", "transcript_id", "model_signature", "full_length_support",
        "evidence_available", "exact_model_matches", "same_intron_chain_models",
        "shared_introns", "shared_exons", "classification",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(model_presence_rows(source_a, index_b, pair_id, "source_a", evidence_a))
        writer.writerows(model_presence_rows(source_b, index_a, pair_id, "source_b", evidence_b))


def write_feature_presence(path: str | Path, pair_id: str, left: dict[str, object], right: dict[str, object]) -> None:
    header = ["pair_id", "feature_type", "feature_signature", "source_a_count", "source_b_count", "classification", "source_a_examples", "source_b_examples"]
    feature_types = ("exon", "intron", "intron_chain", "terminal_exon")
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for feature_type in feature_types:
            keys = sorted(set(left["features"][feature_type]) | set(right["features"][feature_type]))
            for key in keys:
                left_count = left["features"][feature_type].get(key, 0)
                right_count = right["features"][feature_type].get(key, 0)
                classification = "shared" if left_count and right_count else ("source_a_only" if left_count else "source_b_only")
                writer.writerow([
                    pair_id, feature_type, key, left_count, right_count, classification,
                    ",".join(left["examples"].get(key, [])), ",".join(right["examples"].get(key, [])),
                ])


def write_terminal_comparison(path: str | Path, pair_id: str, left: dict[str, object], right: dict[str, object]) -> None:
    header = ["pair_id", "side", "terminal_exon_signature", "source_a_count", "source_b_count", "classification", "source_a_examples", "source_b_examples"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for side in ("5p", "3p"):
            keys = sorted(key for key in set(left["features"]["terminal_exon"]) | set(right["features"]["terminal_exon"]) if key.startswith(f"{side}|"))
            for key in keys:
                left_count = left["features"]["terminal_exon"].get(key, 0)
                right_count = right["features"]["terminal_exon"].get(key, 0)
                classification = "shared" if left_count and right_count else ("source_a_only" if left_count else "source_b_only")
                writer.writerow([
                    pair_id, side, key, left_count, right_count, classification,
                    ",".join(left["examples"].get(key, [])), ",".join(right["examples"].get(key, [])),
                ])


def exon_length(row: dict[str, str]) -> int:
    return sum(end - start + 1 for start, end in intervals(row))


def overlap_length(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    return sum(max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
               for a_start, a_end in left for b_start, b_end in right)


def reciprocal_overlap(left: dict[str, str], right: dict[str, str]) -> tuple[int, float, float]:
    overlap = overlap_length(intervals(left), intervals(right))
    left_len, right_len = exon_length(left), exon_length(right)
    return overlap, overlap / left_len if left_len else 0.0, overlap / right_len if right_len else 0.0


def related(left: dict[str, str], right: dict[str, str], minimum_bp: int = 10,
            reciprocal: float = 0.10) -> bool:
    if left.get("assembly") != right.get("assembly") or left.get("contig") != right.get("contig"):
        return False
    if left.get("strand", ".") != right.get("strand", "."):
        return False
    if junctions(left) & junctions(right):
        return True
    overlap, left_fraction, right_fraction = reciprocal_overlap(left, right)
    return overlap >= minimum_bp and left_fraction >= reciprocal and right_fraction >= reciprocal


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reference_candidates(row: dict[str, str], reference_rows: list[dict[str, str]], minimum_bp: int,
                         reciprocal: float,
                         reference_index: dict[tuple[str, str, str], tuple[list[int], list[tuple[int, int, int]]]]) -> list[tuple[tuple[int, int, float, float, str], dict[str, str]]]:
    candidates = []
    for ref_index in indexed_candidates(row, reference_index):
        ref = reference_rows[ref_index]
        if not valid_gene_id(ref.get("original_gene_id")) or not related(row, ref, minimum_bp, reciprocal):
            continue
        overlap, left_fraction, right_fraction = reciprocal_overlap(row, ref)
        shared_junctions = len(junctions(row) & junctions(ref))
        score = (int(shared_junctions > 0), shared_junctions, overlap, min(left_fraction, right_fraction), ref["original_gene_id"])
        candidates.append((score, ref))
    return sorted(candidates, key=lambda item: item[0], reverse=True)


def assign_loci(source_a: list[dict[str, str]], source_b: list[dict[str, str]], assembly: str,
                reference_rows: list[dict[str, str]] | None = None, minimum_bp: int = 10,
                reciprocal: float = 0.10) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Assign source, reference, then union-derived loci deterministically."""
    source_names = (source_a[0].get("source", "source_a") if source_a else "source_a",
                    source_b[0].get("source", "source_b") if source_b else "source_b")
    all_rows = [empty_assignment(row) for row in source_a + source_b]
    reference_rows = reference_rows or []
    reference_index = build_interval_index(reference_rows)

    for row in all_rows:
        row["assembly"] = row.get("assembly") or assembly
        if valid_gene_id(row.get("original_gene_id")):
            row["gene_assignment_type"] = "biological_source_gene"
            row["locus_id"] = f"source:{row['source']}:{row['original_gene_id']}"
            row["locus_assignment_method"] = "source_gene_id"
        else:
            row["gene_assignment_type"] = "unassigned"
            row["locus_id"] = ""
            row["locus_assignment_method"] = ""

    for row in all_rows:
        if row["locus_id"] or not reference_rows:
            continue
        matches = reference_candidates(row, reference_rows, minimum_bp, reciprocal, reference_index)
        if matches:
            gene_id = matches[0][1]["original_gene_id"]
            row["gene_assignment_type"] = "reference_gene"
            row["locus_id"] = f"reference:{gene_id}"
            row["locus_assignment_method"] = "reference_splice_or_reciprocal_exon_overlap"

    synthetic = [row for row in all_rows if not row["locus_id"] and row.get("model_status") in {"valid", "no_exon_records"}]
    # Build connected components with a sweep-line over overlapping spans.
    # The previous implementation repeatedly scanned all prior interval
    # entries for each model and selected each seed with min(unseen); both
    # become quadratic for millions of independent models. Shared junctions
    # are caught by the same sweep because transcripts sharing a junction also
    # have overlapping transcript spans.
    synthetic_order = sorted(range(len(synthetic)), key=lambda index: (
            synthetic[index].get("assembly", assembly), synthetic[index].get("contig", ""),
            synthetic[index].get("strand", "."), safe_int(synthetic[index].get("transcript_start")),
            synthetic[index].get("transcript_id", "")))
    parent = list(range(len(synthetic)))
    component_size = [1] * len(synthetic)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root == right_root:
            return
        if component_size[left_root] < component_size[right_root]:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        component_size[left_root] += component_size[right_root]

    active: set[int] = set()
    expiry: list[tuple[int, int]] = []
    active_key: tuple[str, str, str] | None = None
    for index in synthetic_order:
        row = synthetic[index]
        span = coordinate_span(row)
        key = (row.get("assembly", assembly), row.get("contig", ""), row.get("strand", "."))
        if key != active_key:
            active.clear()
            expiry.clear()
            active_key = key
        if span is None:
            continue
        start, end = span
        while expiry and expiry[0][0] < start:
            _expired_end, expired = heapq.heappop(expiry)
            active.discard(expired)
        for candidate in active:
            if related(row, synthetic[candidate], minimum_bp, reciprocal):
                union(index, candidate)
        active.add(index)
        heapq.heappush(expiry, (end, index))

    # The first row encountered for a component is its minimum deterministic
    # sort key, matching the component ordering used by the old implementation.
    root_loci: dict[int, str] = {}
    locus_number = 0
    for seed in synthetic_order:
        root = find(seed)
        locus = root_loci.get(root)
        if locus is None:
            locus_number += 1
            first = synthetic[seed]
            locus = f"locus:{assembly}:{first.get('contig', '')}:{first.get('strand', '.')}:" + f"{locus_number:06d}"
            root_loci[root] = locus
        synthetic[seed]["gene_assignment_type"] = "synthetic_locus"
        synthetic[seed]["locus_id"] = locus
        synthetic[seed]["locus_assignment_method"] = "de_novo_union_junction_or_reciprocal_exon_overlap"

    for row in all_rows:
        if not row["locus_id"]:
            row["gene_assignment_type"] = "unassigned"
            row["locus_assignment_method"] = "unassigned_missing_coordinates"

    all_rows.sort(key=lambda row: (row.get("source", ""), row.get("contig", ""), safe_int(row.get("transcript_start")), row.get("transcript_id", "")))
    return ([row for row in all_rows if row.get("source") == source_names[0]],
            [row for row in all_rows if row.get("source") == source_names[1]])


def sha256_if_file(path: str) -> str:
    return sha256(path) if path and Path(path).is_file() and Path(path).stat().st_size else ""


def assignment_summary(rows: list[dict[str, str]], pair_id: str, assembly: str, reference_path: str = "") -> list[dict[str, str]]:
    counts = Counter((row.get("source", ""), row.get("gene_assignment_type", "unassigned"), row.get("locus_assignment_method", "")) for row in rows)
    total_loci = len({row.get("locus_id") for row in rows if row.get("locus_id")})
    reference_checksum = sha256_if_file(reference_path)
    output = []
    for (source, assignment_type, method), count in sorted(counts.items()):
        output.append({
            "pair_id": pair_id, "assembly": assembly, "source": source,
            "gene_assignment_type": assignment_type, "locus_assignment_method": method,
            "n_models": str(count), "n_loci": str(len({row.get("locus_id") for row in rows if row.get("source") == source and row.get("gene_assignment_type") == assignment_type and row.get("locus_id")})),
            "total_pair_loci": str(total_loci), "reference_annotation": reference_path,
            "reference_sha256": reference_checksum, "minimum_exonic_overlap_bp": "10",
            "minimum_reciprocal_exon_overlap": "0.10",
        })
    return output


def cmd_validate(args: argparse.Namespace) -> None:
    errors: list[str] = []
    if args.source_a_assembly and args.source_b_assembly and args.source_a_assembly != args.source_b_assembly:
        errors.append(f"source assemblies differ: {args.source_a_assembly} != {args.source_b_assembly}")
    for label, value in (("source_a_path", args.source_a_path), ("source_b_path", args.source_b_path)):
        if not value or value.lower() in {"na", "null"}:
            errors.append(f"{label} is missing")
    for label, fmt in (("source_a_format", args.source_a_format), ("source_b_format", args.source_b_format)):
        if fmt.lower() not in VALID_FORMATS:
            errors.append(f"{label}={fmt} is not one of gff3, gtf, bed12, auto")
    for label, value in (("source_a_path", args.source_a_path), ("source_b_path", args.source_b_path), ("reference_annotation", args.reference_annotation)):
        if value and value.lower() not in {"na", "null"} and not re.match(r"(?i)^(https?|ftp)://", value) and not Path(value).exists():
            errors.append(f"{label} does not exist: {value}")
    if errors:
        print("pair manifest validation failed: " + "; ".join(errors), file=sys.stderr)
        raise SystemExit(2)
    write_tsv(args.output, ["pair_id", "assembly", "source_a", "source_b", "format_a", "format_b", "reference_annotation"], [[args.pair_id, args.assembly, args.source_a, args.source_b, args.source_a_format.lower(), args.source_b_format.lower(), args.reference_annotation or ""]])


def cmd_normalise(args: argparse.Namespace) -> None:
    rows, warnings = parse_models(args.input, args.source, args.pair_id, args.assembly, args.format)
    write_tsv(args.output, MODEL_HEADER, rows)
    if warnings:
        Path(args.output + ".warnings").write_text("\n".join(warnings) + "\n", encoding="utf-8")


def cmd_stage(args: argparse.Namespace) -> None:
    for value, destination in (
        (args.source_a, args.output_source_a),
        (args.source_b, args.output_source_b),
        (args.assembly_report, args.output_assembly_report),
        (args.evidence_a, args.output_evidence_a),
        (args.evidence_b, args.output_evidence_b),
        (args.reference, args.output_reference),
    ):
        stage_input(value, destination, args.cache_dir)


def cmd_normalise_pair(args: argparse.Namespace) -> None:
    source_a, warnings_a = parse_models(args.source_a, "source_a", args.pair_id, args.assembly, args.format_a)
    source_b, warnings_b = parse_models(args.source_b, "source_b", args.pair_id, args.assembly, args.format_b)
    write_tsv(args.output_a, MODEL_HEADER, source_a)
    write_tsv(args.output_b, MODEL_HEADER, source_b)
    write_tsv(
        args.audit,
        ["pair_id", "source", "n_models", "n_warnings", "status"],
        [
            [args.pair_id, "source_a", len(source_a), len(warnings_a), "warning" if warnings_a else "ok"],
            [args.pair_id, "source_b", len(source_b), len(warnings_b), "warning" if warnings_b else "ok"],
        ],
    )


def load_evidence(path: str | Path) -> dict[str, dict[str, str]]:
    return {
        row.get("transcript_id", ""): row
        for row in iter_tsv(path)
        if row.get("transcript_id")
    }


def cmd_compare_features(args: argparse.Namespace) -> None:
    source_a = build_feature_index(args.source_a)
    source_b = build_feature_index(args.source_b)
    evidence_a = load_evidence(args.evidence_a)
    evidence_b = load_evidence(args.evidence_b)
    write_model_presence(args.model_presence, args.pair_id, args.source_a, args.source_b, source_a, source_b, evidence_a, evidence_b)
    write_feature_presence(args.feature_presence, args.pair_id, source_a, source_b)
    write_terminal_comparison(args.terminal_exons, args.pair_id, source_a, source_b)

    classifications_a = Counter()
    classifications_b = Counter()
    for row in iter_tsv(args.source_a):
        classifications_a[classify_model_difference(feature_signature(row), source_b)[0]] += 1
    for row in iter_tsv(args.source_b):
        classifications_b[classify_model_difference(feature_signature(row), source_a)[0]] += 1
    summary_rows = []
    for source, index, classifications in (("source_a", source_a, classifications_a), ("source_b", source_b, classifications_b)):
        for classification in sorted(classifications):
            summary_rows.append([args.pair_id, source, "model", classification, classifications[classification]])
        for feature_type in ("exon", "intron", "intron_chain", "terminal_exon"):
            summary_rows.append([args.pair_id, source, feature_type, "distinct_features", len(index["features"][feature_type])])
    write_tsv(args.summary, ["pair_id", "source", "feature_type", "classification", "count"], summary_rows)


def cmd_assign(args: argparse.Namespace) -> None:
    source_a = read_tsv(args.source_a)
    source_b = read_tsv(args.source_b)
    reference_rows: list[dict[str, str]] = []
    if args.reference and Path(args.reference).is_file() and Path(args.reference).stat().st_size:
        reference_rows, _warnings = parse_models(args.reference, "reference", args.pair_id, args.assembly, args.reference_format)
    assigned_a, assigned_b = assign_loci(source_a, source_b, args.assembly, reference_rows, args.minimum_overlap_bp, args.reciprocal_overlap)
    write_tsv(args.output_a, MODEL_HEADER, assigned_a)
    write_tsv(args.output_b, MODEL_HEADER, assigned_b)
    write_tsv(args.summary, ["pair_id", "assembly", "source", "gene_assignment_type", "locus_assignment_method", "n_models", "n_loci", "total_pair_loci", "reference_annotation", "reference_sha256", "minimum_exonic_overlap_bp", "minimum_reciprocal_exon_overlap"], assignment_summary(assigned_a + assigned_b, args.pair_id, args.assembly, args.reference or ""))


def cmd_audit_loci(args: argparse.Namespace) -> None:
    rows = read_tsv(args.source_a) + read_tsv(args.source_b)
    audit_header = ["pair_id", "source", "assembly", "locus_id", "gene_assignment_type", "locus_assignment_method", "n_models", "transcript_ids"]
    grouped: defaultdict[tuple[str, str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row.get("pair_id", args.pair_id), row.get("source", ""), row.get("assembly", args.assembly), row.get("locus_id", ""), row.get("gene_assignment_type", "unassigned"), row.get("locus_assignment_method", ""))].append(row)
    audit_rows = [{"pair_id": key[0], "source": key[1], "assembly": key[2], "locus_id": key[3], "gene_assignment_type": key[4], "locus_assignment_method": key[5], "n_models": str(len(value)), "transcript_ids": ",".join(sorted(row.get("transcript_id", row.get("original_transcript_id", "")) for row in value))} for key, value in sorted(grouped.items())]
    write_tsv(args.audit, audit_header, audit_rows)
    write_tsv(args.summary, ["pair_id", "assembly", "source", "gene_assignment_type", "locus_assignment_method", "n_models", "n_loci", "total_pair_loci", "reference_annotation", "reference_sha256", "minimum_exonic_overlap_bp", "minimum_reciprocal_exon_overlap"], assignment_summary(rows, args.pair_id, args.assembly, args.reference or ""))


def model_group_key(row: dict[str, str]) -> str:
    if row.get("locus_id"):
        return row["locus_id"]
    if row.get("original_gene_id"):
        return f"unassigned_gene:{row['original_gene_id']}"
    return f"unassigned_model:{row.get('contig', '')}:{row.get('strand', '.')}:" + row.get("transcript_id", row.get("original_transcript_id", ""))


def load_model_groups(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    groups: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[model_group_key(row)].append(row)
    return groups


def comparison_candidates(row: dict[str, str], groups: dict[str, list[dict[str, str]]], all_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    candidates = list(groups.get(model_group_key(row), []))
    if candidates:
        return candidates
    candidates = [candidate for candidate in all_rows if related(row, candidate)]
    return candidates or [candidate for candidate in all_rows if row.get("original_gene_id") and row.get("original_gene_id") == candidate.get("original_gene_id")]


def cmd_intron(args: argparse.Namespace) -> None:
    a, b = read_tsv(args.source_a), read_tsv(args.source_b)
    b_by_locus = load_model_groups(b)
    rows = []
    for tx in a:
        candidates = comparison_candidates(tx, b_by_locus, b)
        best = max(candidates, key=lambda x: (len(junctions(tx) & junctions(x)), reciprocal_overlap(tx, x)[0], x.get("transcript_id", "")), default=None)
        left, right = junctions(tx), junctions(best) if best else set()
        union = len(left | right)
        rows.append([args.pair_id, tx.get("locus_id", ""), tx.get("transcript_id", tx.get("original_transcript_id", "")), best.get("locus_id", "") if best else "", best.get("transcript_id", best.get("original_transcript_id", "")) if best else "", len(left), len(right), len(left & right), f"{len(left & right) / union:.6f}" if union else "1.000000", "exact" if left == right else ("partial" if left & right else "disagree")])
    write_tsv(args.output, ["pair_id", "locus_id", "source_a_transcript_id", "source_b_locus_id", "source_b_transcript_id", "source_a_introns", "source_b_introns", "shared_introns", "jaccard", "classification"], rows)


def cmd_ends(args: argparse.Namespace) -> None:
    a, b = read_tsv(args.source_a), read_tsv(args.source_b)
    b_by_locus = load_model_groups(b)
    rows = []
    for tx in a:
        candidates = comparison_candidates(tx, b_by_locus, b)
        for side in ("five_prime_end", "three_prime_end"):
            if not candidates:
                rows.append([args.pair_id, tx.get("locus_id", ""), side, tx.get("transcript_id", tx.get("original_transcript_id", "")), "", "", "", "missing_partner"])
                continue
            partner = min(candidates, key=lambda x: (abs(safe_int(tx.get(side), 10**18) - safe_int(x.get(side), 10**18)), x.get("transcript_id", "")))
            distance = abs(safe_int(tx.get(side), 10**18) - safe_int(partner.get(side), 10**18)) if tx.get(side, "").isdigit() and partner.get(side, "").isdigit() else ""
            rows.append([args.pair_id, tx.get("locus_id", ""), side, tx.get("transcript_id", tx.get("original_transcript_id", "")), partner.get("locus_id", ""), partner.get("transcript_id", partner.get("original_transcript_id", "")), distance, "within_threshold" if isinstance(distance, int) and distance <= args.threshold else "outside_threshold"])
    write_tsv(args.output, ["pair_id", "locus_id", "end", "source_a_transcript_id", "source_b_locus_id", "source_b_transcript_id", "distance_bp", "classification"], rows)


def cmd_utr(args: argparse.Namespace) -> None:
    a, b = read_tsv(args.source_a), read_tsv(args.source_b)
    b_by_locus = load_model_groups(b)
    rows = []
    for tx in a:
        partners = comparison_candidates(tx, b_by_locus, b)
        partner = partners[0] if partners else None
        available_a, available_b = bool(tx.get("utr_coordinates")), bool(partner and partner.get("utr_coordinates"))
        same = available_a and available_b and tx["utr_coordinates"] == partner["utr_coordinates"]
        rows.append([args.pair_id, tx.get("locus_id", ""), tx.get("transcript_id", tx.get("original_transcript_id", "")), partner.get("transcript_id", partner.get("original_transcript_id", "")) if partner else "", int(available_a), int(available_b), int(same), "missing" if not (available_a and available_b) else ("concordant" if same else "discordant")])
    write_tsv(args.output, ["pair_id", "locus_id", "source_a_transcript_id", "source_b_transcript_id", "source_a_utr_available", "source_b_utr_available", "exact_match", "classification"], rows)


def cmd_summary(args: argparse.Namespace) -> None:
    chain_rows, end_rows, utr_rows = read_tsv(args.intron), read_tsv(args.ends), read_tsv(args.utr)
    summary = [[args.pair_id, len(chain_rows), sum(r.get("classification") == "exact" for r in chain_rows), sum(r.get("classification") == "partial" for r in chain_rows), len(end_rows), sum(r.get("classification") == "within_threshold" for r in end_rows), len(utr_rows), sum(r.get("classification") == "concordant" for r in utr_rows)]]
    write_tsv(args.output, ["pair_id", "n_intron_comparisons", "n_exact_intron_chains", "n_partial_intron_chains", "n_end_comparisons", "n_ends_within_threshold", "n_utr_comparisons", "n_exact_utr_matches"], summary)
    disagreements = [[args.pair_id, r.get("locus_id", ""), r.get("source_a_transcript_id", ""), r.get("source_b_transcript_id", ""), "intron_chain", r.get("classification", "")] for r in chain_rows if r.get("classification") == "disagree"]
    disagreements += [[args.pair_id, r.get("locus_id", ""), r.get("source_a_transcript_id", ""), r.get("source_b_transcript_id", ""), "transcript_end", r.get("classification", "")] for r in end_rows if r.get("classification") == "outside_threshold"]
    write_tsv(args.disagreements, ["pair_id", "locus_id", "source_a_transcript_id", "source_b_transcript_id", "feature", "classification"], disagreements)


def cmd_evidence(args: argparse.Namespace) -> None:
    evidence = {row.get("transcript_id"): row for row in read_tsv(args.evidence)}
    rows = []
    for row in read_tsv(args.models):
        tx = row.get("transcript_id", row.get("original_transcript_id", ""))
        ev = evidence.get(tx, {})
        rows.append([row.get("pair_id", ""), row.get("source", ""), tx, row.get("locus_id", ""), ev.get("read_count", ""), ev.get("junction_support", ""), ev.get("full_length_support", ""), ev.get("source_tool", "")])
    write_tsv(args.output, ["pair_id", "source", "transcript_id", "locus_id", "read_count", "junction_support", "full_length_support", "source_tool"], rows)


def cmd_audit(args: argparse.Namespace) -> None:
    rows = read_tsv(args.models)
    malformed = sum(1 for row in rows if row.get("model_status") not in {"valid", "no_exon_records"})
    write_tsv(args.output, ["pair_id", "source", "source_format", "n_models", "n_malformed", "n_unique_genes", "n_unique_loci", "conversion_loss", "status"], [[args.pair_id, args.source, args.source_format, len(rows), malformed, len({r.get("original_gene_id") for r in rows if r.get("original_gene_id")}), len({r.get("locus_id") for r in rows if r.get("locus_id")}), args.conversion_loss, "ok" if not malformed else "warning"]])
    write_tsv(args.id_map, ["pair_id", "source", "original_gene_id", "original_transcript_id", "canonical_transcript_id", "gene_assignment_type", "locus_id"], [[r.get("pair_id", args.pair_id), r.get("source", args.source), r.get("original_gene_id", ""), r.get("original_transcript_id", r.get("transcript_id", "")), f"{args.source}:{r.get('transcript_id', r.get('original_transcript_id', ''))}", r.get("gene_assignment_type", "unassigned"), r.get("locus_id", "")] for r in rows])


def cmd_counts(args: argparse.Namespace) -> None:
    rows = read_tsv(args.source_a) + read_tsv(args.source_b)
    output = []
    for source in sorted({row.get("source", "") for row in rows}):
        source_rows = [row for row in rows if row.get("source") == source]
        output.append([args.pair_id, source, "source", len(source_rows), len({row.get("transcript_id", row.get("original_transcript_id", "")) for row in source_rows}), sum(int(row.get("exon_count") or 0) for row in source_rows), len({row.get("locus_id") for row in source_rows if row.get("locus_id")}), "enabled"])
    output.append([args.pair_id, "pair", "locus", len({row.get("locus_id") for row in rows if row.get("locus_id")}), len({row.get("transcript_id", row.get("original_transcript_id", "")) for row in rows}), sum(int(row.get("exon_count") or 0) for row in rows), len({row.get("locus_id") for row in rows if row.get("locus_id")}), "enabled"])
    write_tsv(args.output, ["pair_id", "scope", "summary_level", "locus_count", "transcript_count", "exon_count", "distinct_locus_count", "status"], output)


def synthetic_attributes(gene_id: str, tx_id: str, locus_id: str, source_tx: str, source: str, assignment_type: str) -> str:
    values = {
        "gene_id": gene_id, "transcript_id": tx_id, "locus_id": locus_id,
        "source_transcript_id": source_tx, "synthetic_gene": "true", "source": source,
        "gene_assignment_type": assignment_type,
    }
    return "; ".join(f"{key}={value}" for key, value in values.items())


def gff3_attributes(values: dict[str, str]) -> str:
    return ";".join(f"{key}={value}" for key, value in values.items())


def cmd_render_synthetic(args: argparse.Namespace) -> None:
    rows = read_tsv(args.input)
    by_locus: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("locus_id") and row.get("transcript_start") and row.get("transcript_end"):
            by_locus[row["locus_id"]].append(row)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for locus_id in sorted(by_locus):
            models = sorted(by_locus[locus_id], key=lambda row: (safe_int(row.get("transcript_start")), row.get("transcript_id", "")))
            start = min(safe_int(row["transcript_start"]) for row in models)
            end = max(safe_int(row["transcript_end"]) for row in models)
            contig, strand = models[0].get("contig", "."), models[0].get("strand", ".")
            gene_id = locus_id
            gene_attrs = {
                "ID": gene_id, "gene_id": gene_id, "transcript_id": gene_id,
                "locus_id": locus_id, "source_transcript_id": "", "synthetic_gene": "true",
                "source": args.source, "gene_assignment_type": "synthetic_locus",
            }
            out.write(f"{contig}\tpairwise_qc\tgene\t{start}\t{end}\t.\t{strand}\t.\t{gff3_attributes(gene_attrs)}\n")
            for row in models:
                source_tx = row.get("original_transcript_id", row.get("transcript_id", ""))
                assignment_type = row.get("gene_assignment_type", "unassigned")
                tx_attrs = {
                    "ID": source_tx, "Parent": gene_id, "gene_id": gene_id,
                    "transcript_id": source_tx, "locus_id": locus_id,
                    "source_transcript_id": source_tx, "synthetic_gene": "true",
                    "source": args.source, "gene_assignment_type": assignment_type,
                }
                out.write(f"{contig}\tpairwise_qc\ttranscript\t{row['transcript_start']}\t{row['transcript_end']}\t.\t{strand}\t.\t{gff3_attributes(tx_attrs)}\n")
                for exon_index, exon in enumerate(intervals(row), 1):
                    exon_attrs = {"ID": f"{source_tx}.exon{exon_index}", "Parent": source_tx, **{key: value for key, value in tx_attrs.items() if key not in {"ID", "Parent"}}}
                    out.write(f"{contig}\tpairwise_qc\texon\t{exon[0]}\t{exon[1]}\t.\t{strand}\t.\t{gff3_attributes(exon_attrs)}\n")


def cmd_bed12gtf(args: argparse.Namespace) -> None:
    with open_text(args.input) as handle, open(args.output, "w", encoding="utf-8") as out:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise SystemExit(f"BED12 line {line_no} has fewer than 12 columns")
            chrom, start0, end, name, score, strand = fields[:6]
            count = int(fields[9])
            sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            offsets = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            if len(sizes) != count or len(offsets) != count:
                raise SystemExit(f"BED12 line {line_no} blockCount does not match block lists")
            tx = name or f"bed12_tx_{line_no}"
            gene = f"synthetic_gene:{tx}"
            tx_attrs = f'gene_id "{gene}"; transcript_id "{tx}"; source_format "BED12"; synthetic_gene "true";'
            out.write(f"{chrom}\tbed12_to_gtf\ttranscript\t{int(start0)+1}\t{end}\t{score}\t{strand or '.'}\t.\t{tx_attrs}\n")
            for size, offset in zip(sizes, offsets):
                start = int(start0) + offset + 1
                out.write(f"{chrom}\tbed12_to_gtf\texon\t{start}\t{start + size - 1}\t{score}\t{strand or '.'}\t.\t{tx_attrs}\n")


def cmd_models2gtf(args: argparse.Namespace) -> None:
    with open(args.input, encoding="utf-8") as handle, open(args.output, "w", encoding="utf-8") as out:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row.get("transcript_start") or not row.get("transcript_end"):
                continue
            tx = row.get("transcript_id", row.get("original_transcript_id", ""))
            gene = row.get("locus_id") or row.get("original_gene_id") or f"unassigned:{tx}"
            tags = f'gene_id "{gene}"; transcript_id "{tx}"; locus_id "{row.get("locus_id", "")}"; source_transcript_id "{row.get("original_transcript_id", tx)}"; synthetic_gene "true"; source "{row.get("source", "")}";'
            out.write(f'{row["contig"]}\tpairwise_qc\ttranscript\t{row["transcript_start"]}\t{row["transcript_end"]}\t.\t{row["strand"]}\t.\t{tags}\n')
            for exon in intervals(row):
                out.write(f'{row["contig"]}\tpairwise_qc\texon\t{exon[0]}\t{exon[1]}\t.\t{row["strand"]}\t.\t{tags}\n')


def cmd_status(args: argparse.Namespace) -> None:
    valid_statuses = {"enabled", "not_applicable", "missing_input", "failed"}
    if any(status not in valid_statuses for _metric, status, _reason in args.entry):
        raise SystemExit("metric status must be enabled, not_applicable, missing_input, or failed")
    write_tsv(args.output, ["pair_id", "metric", "status", "reason"], [[args.pair_id, metric, status, reason] for metric, status, reason in args.entry])


def cmd_report(args: argparse.Namespace) -> None:
    summary = read_tsv(args.summary)
    row = summary[0] if summary else {}
    Path(args.output).write_text("# Pairwise annotation QC report\n\n" + f"- Pair: `{args.pair_id}`\n- Intron-chain comparisons: {row.get('n_intron_comparisons', 0)}\n- Exact intron chains: {row.get('n_exact_intron_chains', 0)}\n- End comparisons: {row.get('n_end_comparisons', 0)}\n- UTR comparisons: {row.get('n_utr_comparisons', 0)}\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    v = sub.add_parser("validate")
    for name in ("pair-id", "assembly", "source-a", "source-b", "source-a-path", "source-b-path", "source-a-format", "source-b-format"):
        v.add_argument(f"--{name}", required=True)
    v.add_argument("--source-a-assembly", default=""); v.add_argument("--source-b-assembly", default="")
    v.add_argument("--reference-annotation", default=""); v.add_argument("--output", required=True); v.set_defaults(func=cmd_validate)
    n = sub.add_parser("normalise")
    n.add_argument("--input", required=True); n.add_argument("--source", required=True); n.add_argument("--pair-id", required=True); n.add_argument("--assembly", required=True); n.add_argument("--format", required=True); n.add_argument("--output", required=True); n.set_defaults(func=cmd_normalise)
    stg = sub.add_parser("stage")
    stg.add_argument("--source-a", default=""); stg.add_argument("--source-b", default="")
    stg.add_argument("--assembly-report", default=""); stg.add_argument("--evidence-a", default=""); stg.add_argument("--evidence-b", default=""); stg.add_argument("--reference", default="")
    stg.add_argument("--cache-dir", required=True)
    stg.add_argument("--output-source-a", required=True); stg.add_argument("--output-source-b", required=True); stg.add_argument("--output-assembly-report", required=True)
    stg.add_argument("--output-evidence-a", required=True); stg.add_argument("--output-evidence-b", required=True); stg.add_argument("--output-reference", required=True); stg.set_defaults(func=cmd_stage)
    np = sub.add_parser("normalise-pair")
    np.add_argument("--source-a", required=True); np.add_argument("--source-b", required=True); np.add_argument("--format-a", required=True); np.add_argument("--format-b", required=True)
    np.add_argument("--pair-id", required=True); np.add_argument("--assembly", required=True); np.add_argument("--output-a", required=True); np.add_argument("--output-b", required=True); np.add_argument("--audit", required=True); np.set_defaults(func=cmd_normalise_pair)
    cf = sub.add_parser("compare-features")
    cf.add_argument("--source-a", required=True); cf.add_argument("--source-b", required=True); cf.add_argument("--evidence-a", default=""); cf.add_argument("--evidence-b", default=""); cf.add_argument("--pair-id", required=True)
    cf.add_argument("--model-presence", required=True); cf.add_argument("--feature-presence", required=True); cf.add_argument("--terminal-exons", required=True); cf.add_argument("--summary", required=True); cf.set_defaults(func=cmd_compare_features)
    a = sub.add_parser("assign-loci")
    a.add_argument("--source-a", required=True); a.add_argument("--source-b", required=True); a.add_argument("--assembly", required=True); a.add_argument("--pair-id", required=True); a.add_argument("--reference", default=""); a.add_argument("--reference-format", default="auto"); a.add_argument("--minimum-overlap-bp", type=int, default=10); a.add_argument("--reciprocal-overlap", type=float, default=0.10); a.add_argument("--output-a", required=True); a.add_argument("--output-b", required=True); a.add_argument("--summary", required=True); a.set_defaults(func=cmd_assign)
    al = sub.add_parser("audit-loci")
    al.add_argument("--source-a", required=True); al.add_argument("--source-b", required=True); al.add_argument("--pair-id", required=True); al.add_argument("--assembly", required=True); al.add_argument("--reference", default=""); al.add_argument("--audit", required=True); al.add_argument("--summary", required=True); al.set_defaults(func=cmd_audit_loci)
    au = sub.add_parser("audit")
    au.add_argument("--models", required=True); au.add_argument("--pair-id", required=True); au.add_argument("--source", required=True); au.add_argument("--source-format", required=True); au.add_argument("--conversion-loss", default="0"); au.add_argument("--output", required=True); au.add_argument("--id-map", required=True); au.set_defaults(func=cmd_audit)
    for name, fn in (("intron", cmd_intron), ("ends", cmd_ends), ("utr", cmd_utr)):
        q = sub.add_parser(name); q.add_argument("--source-a", required=True); q.add_argument("--source-b", required=True); q.add_argument("--pair-id", required=True); q.add_argument("--output", required=True); q.set_defaults(func=fn)
    sub.choices["ends"].add_argument("--threshold", type=int, default=50)
    s = sub.add_parser("summary"); s.add_argument("--pair-id", required=True); s.add_argument("--intron", required=True); s.add_argument("--ends", required=True); s.add_argument("--utr", required=True); s.add_argument("--output", required=True); s.add_argument("--disagreements", required=True); s.set_defaults(func=cmd_summary)
    e = sub.add_parser("evidence"); e.add_argument("--models", required=True); e.add_argument("--evidence", required=True); e.add_argument("--output", required=True); e.set_defaults(func=cmd_evidence)
    c = sub.add_parser("counts"); c.add_argument("--source-a", required=True); c.add_argument("--source-b", required=True); c.add_argument("--pair-id", required=True); c.add_argument("--output", required=True); c.set_defaults(func=cmd_counts)
    r = sub.add_parser("render-synthetic"); r.add_argument("--input", required=True); r.add_argument("--source", required=True); r.add_argument("--output", required=True); r.set_defaults(func=cmd_render_synthetic)
    b = sub.add_parser("bed12gtf"); b.add_argument("--input", required=True); b.add_argument("--output", required=True); b.set_defaults(func=cmd_bed12gtf)
    m = sub.add_parser("models2gtf"); m.add_argument("--input", required=True); m.add_argument("--output", required=True); m.set_defaults(func=cmd_models2gtf)
    st = sub.add_parser("status"); st.add_argument("--pair-id", required=True); st.add_argument("--output", required=True); st.add_argument("--entry", nargs=3, action="append", required=True); st.set_defaults(func=cmd_status)
    rp = sub.add_parser("report"); rp.add_argument("--pair-id", required=True); rp.add_argument("--summary", required=True); rp.add_argument("--output", required=True); rp.set_defaults(func=cmd_report)
    return p


if __name__ == "__main__":
    parsed = parser().parse_args()
    parsed.func(parsed)
