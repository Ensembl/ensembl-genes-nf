#!/usr/bin/env python3
"""
Build a Translon call database from mixed translon caller outputs.

The output model is deliberately database-shaped rather than BED-shaped:
`translons` stores one row per called translon, `translon_blocks` stores
child genomic blocks, `parser_manifest` records how each input was interpreted,
and optional CDS recall tables compare calls with CDS intervals from a GTF.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, TextIO

import pandas as pd


STOP_CODONS = {"TAA", "TAG", "TGA"}
NEAR_COGNATE_STARTS = {"CTG", "GTG", "TTG", "ACG", "AGG", "ATA", "ATT", "ATC", "AAG"}
SUPPORTED_SUFFIXES = {".bed", ".bed12", ".gtf", ".gff", ".gff3", ".csv", ".tsv", ".txt"}
PRICE_SAMPLE_MAP = {
    "Fibo": "Ribo_Fib_pooled",
    "Fibo_0": "SRR15513179_1",
    "Fibo_1": "SRR15513180_1",
    "Fibo_2": "SRR15513181_1",
    "Fibo_3": "SRR15513182_1",
    "Fibo_4": "Fib_24_45m",
    "Fibo_5": "Fib_24_bsl",
    "Fibo_6": "Fib_27_45m",
    "Fibo_7": "Fib_27_bsl",
    "Fibo_8": "Fib_41_45m",
    "Fibo_9": "Fib_41_bsl",
    "Endo": "Ribo_EC_pooled",
    "Endo_0": "SRR15513197",
    "Endo_1": "SRR15513198_GENELAB1026",
    "Endo_2": "SRR15513199",
    "Endo_3": "SRR15513200",
    "Endo_4": "SRR15513201",
    "Endo_5": "SRR15513202_GENELAB1143",
    "Pancreas": "Ribo_Pancreas_pooled",
    "Pancreas_0": "SRR11005875_to_79",
    "Pancreas_1": "SRR11005880_to_84",
    "Pancreas_2": "SRR11005885_to_89",
    "Pancreas_3": "SRR11005890_to_94",
    "Pancreas_4": "SRR11005895_to_99",
    "Pancreas_5": "SRR11005900_to_04",
    "pancreas_mymapping": "Ribo_Pancreas_pooled_fastq",
    "pancreas_mymapping_0": "SRR11005875_to_79_fastq",
    "pancreas_mymapping_1": "SRR11005880_to_84_fastq",
    "pancreas_mymapping_2": "SRR11005885_to_89_fastq",
    "pancreas_mymapping_3": "SRR11005890_to_94_fastq",
    "pancreas_mymapping_4": "SRR11005895_to_99_fastq",
    "pancreas_mymapping_5": "SRR11005900_to_04_fastq",
}
CONVERSION_RULES = {
    "orfquant_gff": "orfquant_bedlike_plus_stop",
    "ribotie_gtf": "ribotie_gff1_plus_stop",
    "iribo_gff": "iribo_bedlike",
    "generic_gff": "gff_1based",
    "bed12": "bed12_native",
    "translonscorer_csv": "csv_interval_native",
}
EMPTY_TABLE_COLUMNS = {
    "reference_cds": [
        "transcript_id",
        "gene_id",
        "gene_name",
        "bed_chrom",
        "bed_start",
        "bed_end",
        "bed_strand",
        "block_count",
        "block_sizes",
        "block_starts",
        "spliced_length_nt",
        "feature_key",
    ],
    "cds_recall_by_tool": [
        "source_tool",
        "reference_cds",
        "exact_cds_recalled",
        "exact_cds_recall_pct",
    ],
}


@dataclass
class Candidate:
    source_tool: str
    parser_name: str
    raw_file: str
    raw_label: str
    sample_id: str
    native_translon_id: str
    chrom: str
    strand: str
    intervals: list[tuple[int, int]]
    score: str = "0"
    transcript_id: str = ""
    gene_id: str = ""
    gene_name: str = ""
    native_feature_type: str = ""
    attributes: dict[str, str] = field(default_factory=dict)


@dataclass
class InputRecord:
    path: Path
    tool_hint: str
    sample_id: str


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")  # type: ignore[return-value]
    return path.open()


def normalise_sample(raw: str) -> str:
    sample = raw.replace("GENELAB-000", "GENELAB").replace("GENELAB_000", "GENELAB")
    sample = re.sub(r"\.(known|novel)$", "", sample)
    sample = PRICE_SAMPLE_MAP.get(sample, sample)
    return re.sub(r"_1$", "", sample)


def chrom_to_ucsc(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"MT", "M"}:
        return "chrM"
    return f"chr{chrom}"


def ucsc_to_seq_region(chrom: str) -> str:
    if chrom == "chrM":
        return "MT"
    return chrom[3:] if chrom.startswith("chr") else chrom


def parse_attrs(attrs: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+)\s+"([^"]+)"', attrs):
        parsed[key] = value
    for part in attrs.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            parsed[key.strip()] = value.strip().strip('"')
        elif " " in part:
            key, value = part.split(None, 1)
            parsed[key.strip()] = value.strip().strip('"')
    return parsed


def split_int_list(value: str) -> list[int]:
    return [int(part) for part in value.rstrip(",").split(",") if part]


def detect_parser(path: Path, tool_hint: str = "") -> tuple[str, str, str]:
    """Return (status, parser_name, detected_tool)."""
    try:
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#") or line.startswith("track"):
                    continue
                fields = line.rstrip("\n").split("\t")
                source = fields[1] if len(fields) > 1 else ""
                attrs = fields[8] if len(fields) >= 9 else ""
                name = fields[3] if len(fields) >= 4 else ""
                lower_name = path.name.lower()

                if len(fields) >= 12:
                    try:
                        int(fields[1]); int(fields[2]); int(fields[9])
                    except ValueError:
                        pass
                    else:
                        tool = "PRICE" if "price" in lower_name else ""
                        tool = tool or (tool_hint if tool_hint and tool_hint not in {"known", "novel"} else "")
                        tool = tool or "BED12"
                        return "ok", "bed12", tool

                if len(fields) >= 9:
                    if source == "RiboTIE" or "ORF_id" in attrs or "ribotie" in lower_name:
                        return "ok", "ribotie_gtf", "RiboTIE"
                    if source == "iRibo" or "iribo" in lower_name:
                        return "ok", "iribo_gff", "iRibo"
                    if source == "ORFQuant" or ("ENST" in attrs and "=" not in attrs) or "orfquant" in lower_name:
                        return "ok", "orfquant_gff", "ORFQuant"
                    return "ok", "generic_gff", tool_hint or source or "unknown"

                if path.suffix.lower() == ".csv" or "," in line:
                    return "ok", "translonscorer_csv", tool_hint or "TranslonScorer"

                if name:
                    return "unsupported", "unknown", tool_hint or "unknown"
                break
    except Exception:
        return "error", "unreadable", tool_hint or "unknown"
    return "unsupported", "empty_or_unknown", tool_hint or "unknown"


def converted_intervals(
    raw_intervals: list[tuple[int, int]],
    strand: str,
    parser_name: str,
) -> tuple[list[tuple[int, int]], str]:
    rule = CONVERSION_RULES.get(parser_name, "unknown")
    intervals = list(raw_intervals)
    if rule == "ribotie_gff1_plus_stop" or rule == "gff_1based":
        intervals = [(start - 1, end) for start, end in intervals]
    if rule in {"orfquant_bedlike_plus_stop", "ribotie_gff1_plus_stop"}:
        if strand == "+":
            idx = max(range(len(intervals)), key=lambda i: intervals[i][1])
            start, end = intervals[idx]
            intervals[idx] = (start, end + 3)
        elif strand == "-":
            idx = min(range(len(intervals)), key=lambda i: intervals[i][0])
            start, end = intervals[idx]
            intervals[idx] = (start - 3, end)
    return sorted(intervals), rule


def discover_inputs(input_root: Path, manifest: Path | None) -> list[InputRecord]:
    if manifest:
        rows: list[InputRecord] = []
        with manifest.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                path = Path(row["path"])
                rows.append(
                    InputRecord(
                        path=path,
                        tool_hint=row.get("tool", "") or path.parent.name,
                        sample_id=row.get("sample_id", "") or normalise_sample(path.stem),
                    )
                )
        return rows

    records = []
    for path in sorted(input_root.rglob("*")):
        suffixes = "".join(path.suffixes[-2:]).lower()
        suffix = path.suffix.lower()
        if path.is_file() and (suffix in SUPPORTED_SUFFIXES or suffixes.endswith(".gtf.gz") or suffixes.endswith(".bed.gz")):
            records.append(InputRecord(path=path, tool_hint=path.parent.name, sample_id=normalise_sample(path.stem)))
    return records


def group_gff(path: Path, parser_name: str, source_tool: str, sample_id: str) -> Iterable[Candidate]:
    grouped: dict[str, list[tuple[str, str, int, int, dict[str, str], str]]] = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            feature = fields[2]
            if feature not in {"CDS", "exon", "translon", "orf"}:
                continue
            if parser_name == "orfquant_gff":
                native_id = fields[8].strip()
            elif parser_name == "ribotie_gtf":
                native_id = attrs.get("ORF_id") or attrs.get("transcript_id") or fields[8].strip()
            else:
                native_id = attrs.get("ID") or attrs.get("Parent") or attrs.get("transcript_id") or fields[8].strip()
            grouped[native_id].append(
                (chrom_to_ucsc(fields[0]), fields[6].strip(), int(fields[3]), int(fields[4]), attrs, fields[5])
            )

    for native_id, rows in grouped.items():
        chroms = {row[0] for row in rows}
        strands = {row[1] for row in rows}
        if len(chroms) != 1 or len(strands) != 1:
            continue
        attrs = rows[0][4]
        intervals, conversion_rule = converted_intervals(
            [(row[2], row[3]) for row in rows],
            next(iter(strands)),
            parser_name,
        )
        yield Candidate(
            source_tool=source_tool,
            parser_name=parser_name,
            raw_file=str(path),
            raw_label=path.stem,
            sample_id=sample_id,
            native_translon_id=native_id,
            chrom=next(iter(chroms)),
            strand=next(iter(strands)),
            intervals=intervals,
            score=rows[0][5] if rows[0][5] != "." else "0",
            transcript_id=attrs.get("transcript_id", "") or parse_enst(native_id),
            gene_id=attrs.get("gene_id", ""),
            gene_name=attrs.get("gene_name", ""),
            native_feature_type=attrs.get("translon_type", attrs.get("type", "")),
            attributes={**attrs, "conversion_rule": conversion_rule},
        )


def parse_enst(value: str) -> str:
    match = re.search(r"(ENST\d+(?:\.\d+)?)", value)
    return match.group(1) if match else ""


def parse_bed12(path: Path, parser_name: str, source_tool: str, sample_id: str) -> Iterable[Candidate]:
    with open_text(path) as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            chrom_start = int(fields[1])
            block_count = int(fields[9])
            sizes = split_int_list(fields[10])
            starts = split_int_list(fields[11])
            if block_count != len(sizes) or block_count != len(starts):
                continue
            intervals = [(chrom_start + rel, chrom_start + rel + size) for rel, size in zip(starts, sizes)]
            yield Candidate(
                source_tool=source_tool,
                parser_name=parser_name,
                raw_file=str(path),
                raw_label=path.stem,
                sample_id=sample_id,
                native_translon_id=fields[3] or f"{path.stem}_{idx}",
                chrom=chrom_to_ucsc(fields[0]),
                strand=fields[5],
                intervals=sorted(intervals),
                score=fields[4],
                transcript_id=parse_enst(fields[3]),
            )


def parse_translonscorer_csv(path: Path, source_tool: str, sample_id: str) -> Iterable[Candidate]:
    with open_text(path) as handle:
        dialect = csv.Sniffer().sniff(handle.read(4096), delimiters=",\t")
        handle.seek(0)
        reader = csv.DictReader(handle, dialect=dialect)
        for idx, row in enumerate(reader, start=1):
            lower = {key.lower(): key for key in row}
            chrom_key = lower.get("chrom") or lower.get("chromosome") or lower.get("seq_region_name")
            start_key = lower.get("start") or lower.get("chrom_start") or lower.get("seq_region_start")
            end_key = lower.get("end") or lower.get("chrom_end") or lower.get("seq_region_end")
            strand_key = lower.get("strand") or lower.get("seq_region_strand")
            if not (chrom_key and start_key and end_key):
                continue
            start = int(float(row[start_key]))
            end = int(float(row[end_key]))
            if start > 0 and "seq_region" in start_key:
                start -= 1
            name = row.get(lower.get("id", ""), "") or row.get(lower.get("orf_id", ""), "") or f"{path.stem}_{idx}"
            strand = row[strand_key] if strand_key else "+"
            if strand == "1":
                strand = "+"
            elif strand == "-1":
                strand = "-"
            yield Candidate(
                source_tool=source_tool,
                parser_name="translonscorer_csv",
                raw_file=str(path),
                raw_label=path.stem,
                sample_id=sample_id,
                native_translon_id=name,
                chrom=chrom_to_ucsc(row[chrom_key]),
                strand=strand,
                intervals=[(start, end)],
                transcript_id=parse_enst(name),
            )


def sequence_from_blocks(genome, chrom: str, intervals: list[tuple[int, int]], strand: str) -> tuple[str, list[str]]:
    if genome is None:
        return "", ["sequence_not_checked"]
    if chrom not in genome:
        return "", ["missing_chrom"]
    try:
        pieces = [str(genome[chrom][start:end]).upper() for start, end in intervals]
    except Exception:
        return "", ["sequence_fetch_failed"]
    seq = "".join(pieces)
    if strand == "-":
        table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        seq = seq.translate(table)[::-1].upper()
    return seq, []


def feature_key(chrom: str, start: int, end: int, strand: str, sizes: str, starts: str) -> str:
    return f"{chrom}:{start}-{end}:{strand}:{sizes}:{starts}"


def start_class(codon: str) -> str:
    if len(codon) != 3:
        return "missing"
    if codon == "ATG":
        return "ATG"
    if codon in NEAR_COGNATE_STARTS:
        return "near_cognate"
    return "other"


def terminal_class(codon: str) -> str:
    if len(codon) != 3:
        return "missing"
    return "stop" if codon in STOP_CODONS else "non_stop"


def build_translon_tables(candidates: Iterable[Candidate], fasta: Path | None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genome = None
    if fasta:
        try:
            from pyfaidx import Fasta

            genome = Fasta(str(fasta), as_raw=True, sequence_always_upper=True)
        except Exception as exc:
            print(f"[WARN] Could not open FASTA with pyfaidx, skipping sequence QC: {fasta}: {exc}", file=sys.stderr)

    translon_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    qc_counter: Counter[tuple[str, str]] = Counter()

    for candidate in candidates:
        intervals = sorted(candidate.intervals)
        if not intervals:
            continue
        bed_start = min(start for start, _ in intervals)
        bed_end = max(end for _, end in intervals)
        sizes = [end - start for start, end in intervals]
        starts = [start - bed_start for start, _ in intervals]
        block_sizes = ",".join(map(str, sizes))
        block_starts = ",".join(map(str, starts))
        seq, reasons = sequence_from_blocks(genome, candidate.chrom, intervals, candidate.strand)
        spliced_len = sum(sizes)
        start_codon = seq[:3] if len(seq) >= 3 else ""
        terminal_codon = seq[-3:] if len(seq) >= 3 else ""
        if spliced_len <= 0:
            reasons.append("zero_spliced_length")
        if spliced_len % 3 != 0:
            reasons.append("length_not_multiple_of_3")
        if genome is not None and terminal_codon not in STOP_CODONS:
            reasons.append("terminal_non_stop")
        qc_status = "pass" if not reasons else "fail"
        qc_counter[(candidate.source_tool, qc_status)] += 1
        translon_id = f"translon_{len(translon_rows):010d}"
        key = feature_key(candidate.chrom, bed_start, bed_end, candidate.strand, block_sizes, block_starts)

        translon_rows.append(
            {
                "translon_id": translon_id,
                "source_tool": candidate.source_tool,
                "parser_name": candidate.parser_name,
                "raw_file": candidate.raw_file,
                "raw_label": candidate.raw_label,
                "sample_id": candidate.sample_id,
                "native_translon_id": candidate.native_translon_id,
                "transcript_id": candidate.transcript_id,
                "gene_id": candidate.gene_id,
                "gene_name": candidate.gene_name,
                "native_feature_type": candidate.native_feature_type,
                "score": candidate.score,
                "conversion_rule": candidate.attributes.get("conversion_rule", CONVERSION_RULES.get(candidate.parser_name, "")),
                "seq_region_name": ucsc_to_seq_region(candidate.chrom),
                "seq_region_start": bed_start + 1,
                "seq_region_end": bed_end,
                "seq_region_strand": 1 if candidate.strand == "+" else -1,
                "bed_chrom": candidate.chrom,
                "bed_start": bed_start,
                "bed_end": bed_end,
                "bed_strand": candidate.strand,
                "block_count": len(intervals),
                "block_sizes": block_sizes,
                "block_starts": block_starts,
                "spliced_length_nt": spliced_len,
                "length_mod3": spliced_len % 3,
                "start_codon": start_codon,
                "terminal_codon": terminal_codon,
                "start_codon_class": start_class(start_codon),
                "terminal_codon_class": terminal_class(terminal_codon),
                "qc_status": qc_status,
                "qc_fail_reasons": ";".join(sorted(set(reasons))),
                "feature_key": key,
            }
        )

        order = list(range(len(intervals)))
        if candidate.strand == "-":
            order = list(reversed(order))
        translation_rank = {idx: rank + 1 for rank, idx in enumerate(order)}
        for genomic_rank, (start, end) in enumerate(intervals, start=1):
            idx = genomic_rank - 1
            block_rows.append(
                {
                    "translon_id": translon_id,
                    "source_tool": candidate.source_tool,
                    "genomic_block_rank": genomic_rank,
                    "translation_block_rank": translation_rank[idx],
                    "seq_region_name": ucsc_to_seq_region(candidate.chrom),
                    "seq_region_start": start + 1,
                    "seq_region_end": end,
                    "seq_region_strand": 1 if candidate.strand == "+" else -1,
                    "bed_chrom": candidate.chrom,
                    "bed_start": start,
                    "bed_end": end,
                    "block_length_nt": end - start,
                }
            )

    if genome is not None:
        genome.close()

    qc = pd.DataFrame(
        [{"source_tool": tool, "qc_status": status, "translons": count} for (tool, status), count in sorted(qc_counter.items())]
    )
    return pd.DataFrame(translon_rows), pd.DataFrame(block_rows), qc


def reference_cds_from_gtf(gtf: Path | None) -> pd.DataFrame:
    if not gtf:
        return pd.DataFrame()
    grouped: dict[str, dict[str, object]] = {}
    with open_text(gtf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = parse_attrs(fields[8])
            tid = attrs.get("transcript_id", "")
            if not tid:
                continue
            record = grouped.setdefault(
                tid,
                {
                    "transcript_id": tid,
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                    "chrom": chrom_to_ucsc(fields[0]),
                    "strand": fields[6],
                    "intervals": [],
                },
            )
            record["intervals"].append((int(fields[3]) - 1, int(fields[4])))  # type: ignore[index,union-attr]

    rows = []
    for record in grouped.values():
        intervals = sorted(record["intervals"])  # type: ignore[arg-type]
        bed_start = min(start for start, _ in intervals)
        bed_end = max(end for _, end in intervals)
        sizes = ",".join(str(end - start) for start, end in intervals)
        starts = ",".join(str(start - bed_start) for start, _ in intervals)
        rows.append(
            {
                "transcript_id": record["transcript_id"],
                "gene_id": record["gene_id"],
                "gene_name": record["gene_name"],
                "bed_chrom": record["chrom"],
                "bed_start": bed_start,
                "bed_end": bed_end,
                "bed_strand": record["strand"],
                "block_count": len(intervals),
                "block_sizes": sizes,
                "block_starts": starts,
                "spliced_length_nt": sum(end - start for start, end in intervals),
                "feature_key": feature_key(record["chrom"], bed_start, bed_end, record["strand"], sizes, starts),
            }
        )
    return pd.DataFrame(rows)


def cds_recall(translons: pd.DataFrame, reference_cds: pd.DataFrame) -> pd.DataFrame:
    if translons.empty or reference_cds.empty:
        return pd.DataFrame()
    rows = []
    reference_keys = set(reference_cds["feature_key"])
    total = len(reference_keys)
    for tool, group in translons.groupby("source_tool", observed=True):
        called = set(group["feature_key"])
        rows.append(
            {
                "source_tool": tool,
                "reference_cds": total,
                "exact_cds_recalled": len(reference_keys & called),
                "exact_cds_recall_pct": 100 * len(reference_keys & called) / total if total else 0.0,
            }
        )
    called_any = set(translons["feature_key"])
    rows.append(
        {
            "source_tool": "ANY",
            "reference_cds": total,
            "exact_cds_recalled": len(reference_keys & called_any),
            "exact_cds_recall_pct": 100 * len(reference_keys & called_any) / total if total else 0.0,
        }
    )
    return pd.DataFrame(rows)


def write_sqlite(path: Path, tables: dict[str, pd.DataFrame]) -> None:
    if path.exists():
        path.unlink()
    with sqlite3.connect(path) as con:
        for name, df in tables.items():
            if len(df.columns) == 0 and name in EMPTY_TABLE_COLUMNS:
                df = pd.DataFrame(columns=EMPTY_TABLE_COLUMNS[name])
            df.to_sql(name, con, index=False, if_exists="replace")
        con.executescript(
            """
            CREATE INDEX IF NOT EXISTS idx_translons_tool_sample ON translons(source_tool, sample_id);
            CREATE INDEX IF NOT EXISTS idx_translons_feature_key ON translons(feature_key);
            CREATE INDEX IF NOT EXISTS idx_translon_blocks_translon ON translon_blocks(translon_id);
            CREATE INDEX IF NOT EXISTS idx_ref_cds_feature_key ON reference_cds(feature_key);
            """
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--gtf", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    candidates: list[Candidate] = []
    manifest_rows: list[dict[str, object]] = []
    for record in discover_inputs(args.input_root, args.manifest):
        status, parser_name, detected_tool = detect_parser(record.path, record.tool_hint)
        before = len(candidates)
        if status == "ok" and parser_name == "bed12":
            candidates.extend(parse_bed12(record.path, parser_name, detected_tool, record.sample_id))
        elif status == "ok" and parser_name.endswith(("gff", "gtf")):
            candidates.extend(group_gff(record.path, parser_name, detected_tool, record.sample_id))
        elif status == "ok" and parser_name == "translonscorer_csv":
            candidates.extend(parse_translonscorer_csv(record.path, detected_tool, record.sample_id))
        manifest_rows.append(
            {
                "path": str(record.path),
                "tool_hint": record.tool_hint,
                "sample_id": record.sample_id,
                "parser_status": status,
                "parser_name": parser_name,
                "detected_tool": detected_tool,
                "translons": len(candidates) - before,
            }
        )

    parser_manifest = pd.DataFrame(manifest_rows)
    translons, blocks, qc = build_translon_tables(candidates, args.fasta)
    reference_cds = reference_cds_from_gtf(args.gtf)
    recall = cds_recall(translons, reference_cds)

    tables = {
        "translons": translons,
        "translon_blocks": blocks,
        "parser_manifest": parser_manifest,
        "qc_summary": qc,
        "reference_cds": reference_cds,
        "cds_recall_by_tool": recall,
    }
    for name, df in tables.items():
        df.to_csv(args.out_dir / f"{name}.tsv.gz", sep="\t", index=False, compression="gzip")
    write_sqlite(args.out_dir / "translons.sqlite", tables)

    summary = {
        "inputs": int(len(parser_manifest)),
        "parsed_inputs": int(parser_manifest["parser_status"].eq("ok").sum()) if not parser_manifest.empty else 0,
        "translons": int(len(translons)),
        "translon_blocks": int(len(blocks)),
        "reference_cds": int(len(reference_cds)),
    }
    (args.out_dir / "translon_db_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
