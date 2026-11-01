#!/usr/bin/env python3
"""Format-neutral helpers for pairwise annotation QC.

The command deliberately keeps the public contract small.  Nextflow modules
invoke one subcommand per task, while all format parsing and output schemas
live here so GFF3, GTF and BED12 follow the same rules.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote


MODEL_HEADER = [
    "pair_id", "source", "original_gene_id", "original_transcript_id",
    "contig", "strand", "transcript_start", "transcript_end",
    "exon_count", "exon_coordinates", "intron_chain", "five_prime_end",
    "three_prime_end", "assembly", "source_format", "model_status",
    "utr_coordinates", "cds_coordinates",
]


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            match = re.match(r"([^\s]+)\s+\"?([^\"]+)\"?", item)
            if not match:
                continue
            key, value = match.groups()
        result[key.strip()] = unquote(value.strip().strip('"'))
    return result


def fmt_coord(start: int, end: int) -> str:
    return f"{start}-{end}"


def chain(exons: list[tuple[int, int]]) -> str:
    ordered = sorted(exons)
    return ";".join(f"{right}-{next_left}" for (_, right), (next_left, _) in zip(ordered, ordered[1:]))


def infer_format(path: str, requested: str | None = None) -> str:
    if requested and requested.lower() not in {"auto", ""}:
        return requested.lower()
    name = str(path).lower()
    if name.endswith((".bed", ".bed12", ".bed.gz", ".bed12.gz")):
        return "bed12"
    if name.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    return "gff3"


def parse_models(path: str, source: str, pair_id: str, assembly: str, requested: str | None = None) -> tuple[list[dict[str, str]], list[str]]:
    fmt = infer_format(path, requested)
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
                    sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x != ""]
                    offsets = [int(x) for x in fields[11].rstrip(",").split(",") if x != ""]
                    if len(sizes) != block_count or len(offsets) != block_count:
                        raise ValueError("blockCount does not match block lists")
                    exons = [(start0 + offset + 1, start0 + offset + size) for offset, size in zip(offsets, sizes)]
                    if not exons or min(exons)[0] < 1 or max(exons)[1] > end:
                        raise ValueError("invalid BED12 block coordinates")
                except (ValueError, IndexError) as exc:
                    warnings.append(f"bed12_line_{line_no}: {exc}")
                    continue
                tx = name or f"bed12_tx_{line_no}"
                gene = tx.split(".", 1)[0] if tx else f"bed12_gene_{line_no}"
                records[tx] = {"gene": gene, "tx": tx, "contig": chrom, "strand": strand or ".", "exons": exons, "utr": [], "cds": [], "status": "valid"}
    else:
        tx_types = {"transcript", "mRNA", "lnc_RNA", "ncRNA", "miRNA", "snRNA", "snoRNA", "pseudogenic_transcript", "antisense_RNA", "tRNA", "rRNA", "scRNA", "guide_RNA", "vault_RNA", "Y_RNA"}
        pending_exons: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
        pending_utr: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
        pending_cds: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
        tx_gene: dict[str, str] = {}
        tx_meta: dict[str, dict] = {}
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
                if feature in tx_types or feature.lower() in {"transcript", "mrna"}:
                    tx = a.get("transcript_id") or a.get("ID") or a.get("transcript")
                    gene = a.get("gene_id") or a.get("Parent") or a.get("gene") or a.get("gene_name")
                    if tx:
                        tx_gene[tx] = gene or f"synthetic_gene:{tx}"
                        tx_meta[tx] = {"gene": tx_gene[tx], "contig": contig, "strand": strand or ".", "start": start_i, "end": end_i, "status": "valid"}
                elif feature == "exon":
                    parents = (a.get("Parent") or a.get("transcript_id") or a.get("transcript") or "").split(",")
                    for tx in filter(None, parents):
                        pending_exons[tx].append((start_i, end_i))
                        tx_gene.setdefault(tx, a.get("gene_id") or a.get("gene") or f"synthetic_gene:{tx}")
                elif feature.lower() in {"utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr"}:
                    for tx in filter(None, (a.get("Parent") or a.get("transcript_id") or "").split(",")):
                        pending_utr[tx].append((start_i, end_i))
                elif feature == "CDS":
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
            if exons:
                start_i, end_i = min(x[0] for x in exons), max(x[1] for x in exons)
            else:
                start_i = meta.get("start", "")
                end_i = meta.get("end", "")
            records[tx] = {"gene": gene, "tx": tx, "contig": meta.get("contig", ""), "strand": meta.get("strand", "."), "exons": exons, "utr": sorted(pending_utr.get(tx, [])), "cds": sorted(pending_cds.get(tx, [])), "status": status, "start": start_i, "end": end_i}

    rows: list[dict[str, str]] = []
    for rec in records.values():
        exons = sorted(rec.get("exons", []))
        start_i = rec.get("start") or (min(x[0] for x in exons) if exons else "")
        end_i = rec.get("end") or (max(x[1] for x in exons) if exons else "")
        if rec.get("strand") == "-":
            five, three = end_i, start_i
        else:
            five, three = start_i, end_i
        rows.append({
            "pair_id": pair_id, "source": source, "original_gene_id": rec["gene"], "original_transcript_id": rec["tx"],
            "contig": rec.get("contig", ""), "strand": rec.get("strand", "."), "transcript_start": str(start_i), "transcript_end": str(end_i),
            "exon_count": str(len(exons)), "exon_coordinates": ";".join(fmt_coord(*x) for x in exons), "intron_chain": chain(exons),
            "five_prime_end": str(five), "three_prime_end": str(three), "assembly": assembly, "source_format": fmt,
            "model_status": rec.get("status", "valid"), "utr_coordinates": ";".join(fmt_coord(*x) for x in rec.get("utr", [])),
            "cds_coordinates": ";".join(fmt_coord(*x) for x in rec.get("cds", [])),
        })
    return sorted(rows, key=lambda row: (row["contig"], row["transcript_start"], row["original_transcript_id"])), warnings


def read_tsv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str, header: list[str], rows: list[dict | list]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([row.get(column, "") if isinstance(row, dict) else value for column, value in ((column, None) for column in header)] if isinstance(row, dict) else row)


def cmd_validate(args: argparse.Namespace) -> None:
    errors: list[str] = []
    if args.source_a_assembly and args.source_b_assembly and args.source_a_assembly != args.source_b_assembly:
        errors.append(f"source assemblies differ: {args.source_a_assembly} != {args.source_b_assembly}")
    for label, value in (("source_a_path", args.source_a_path), ("source_b_path", args.source_b_path)):
        if not value or value.lower() in {"na", "null"}:
            errors.append(f"{label} is missing")
    for label, fmt in (("source_a_format", args.source_a_format), ("source_b_format", args.source_b_format)):
        if fmt.lower() not in {"gff3", "gtf", "bed12", "auto"}:
            errors.append(f"{label}={fmt} is not one of gff3, gtf, bed12")
    if args.source_a_path and not re.match(r"(?i)^(https?|ftp)://", args.source_a_path) and not Path(args.source_a_path).exists():
        errors.append(f"source_a_path does not exist: {args.source_a_path}")
    if args.source_b_path and not re.match(r"(?i)^(https?|ftp)://", args.source_b_path) and not Path(args.source_b_path).exists():
        errors.append(f"source_b_path does not exist: {args.source_b_path}")
    if errors:
        print("pair manifest validation failed: " + "; ".join(errors), file=sys.stderr)
        raise SystemExit(2)
    write_tsv(args.output, ["pair_id", "assembly", "source_a", "source_b", "format_a", "format_b"], [[args.pair_id, args.assembly, args.source_a, args.source_b, args.source_a_format.lower(), args.source_b_format.lower()]])


def cmd_audit(args: argparse.Namespace) -> None:
    rows = read_tsv(args.models)
    statuses = Counter(row.get("model_status", "") for row in rows)
    malformed = sum(1 for row in rows if row.get("model_status") not in {"valid", "no_exon_records"})
    write_tsv(args.output, ["pair_id", "source", "source_format", "n_models", "n_malformed", "n_unique_genes", "conversion_loss", "status"], [[args.pair_id, args.source, args.source_format, len(rows), malformed, len({r.get('original_gene_id') for r in rows}), args.conversion_loss, "ok" if not malformed else "warning"]])
    write_tsv(args.id_map, ["pair_id", "source", "original_gene_id", "original_transcript_id", "canonical_transcript_id"], [[r["pair_id"], r["source"], r["original_gene_id"], r["original_transcript_id"], f"{args.source}:{r['original_transcript_id']}"] for r in rows])


def load_model_groups(path: str) -> dict[str, list[dict[str, str]]]:
    groups: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(path):
        groups[row.get("original_gene_id", "")].append(row)
    return groups


def cmd_intron(args: argparse.Namespace) -> None:
    a, b = read_tsv(args.source_a), read_tsv(args.source_b)
    b_by_gene = load_model_groups(args.source_b)
    rows = []
    for tx in a:
        candidates = b_by_gene.get(tx["original_gene_id"], [])
        best = max(
            candidates,
            key=lambda x: len(
                set(filter(None, tx["intron_chain"].split(";")))
                & set(filter(None, x["intron_chain"].split(";")))
            ),
            default=None,
        )
        left = set(filter(None, tx["intron_chain"].split(";")))
        right = set(filter(None, best["intron_chain"].split(";"))) if best else set()
        union = len(left | right)
        rows.append([args.pair_id, tx["original_transcript_id"], best["original_transcript_id"] if best else "", len(left), len(right), len(left & right), f"{len(left & right) / union:.6f}" if union else "1.000000", "exact" if left == right else ("partial" if left & right else "disagree")])
    write_tsv(args.output, ["pair_id", "source_a_transcript_id", "source_b_transcript_id", "source_a_introns", "source_b_introns", "shared_introns", "jaccard", "classification"], rows)


def cmd_ends(args: argparse.Namespace) -> None:
    a, b = read_tsv(args.source_a), read_tsv(args.source_b)
    b_by_gene = load_model_groups(args.source_b)
    rows = []
    for tx in a:
        best = b_by_gene.get(tx["original_gene_id"], [])
        for side in ("five_prime_end", "three_prime_end"):
            if not best:
                rows.append([args.pair_id, side, tx["original_transcript_id"], "", "", "missing_partner"])
                continue
            partner = min(best, key=lambda x: abs(int(tx[side]) - int(x[side])) if tx[side].isdigit() and x[side].isdigit() else 10**18)
            distance = abs(int(tx[side]) - int(partner[side])) if tx[side].isdigit() and partner[side].isdigit() else ""
            rows.append([args.pair_id, side, tx["original_transcript_id"], partner["original_transcript_id"], distance, "within_threshold" if isinstance(distance, int) and distance <= args.threshold else "outside_threshold"])
    write_tsv(args.output, ["pair_id", "end", "source_a_transcript_id", "source_b_transcript_id", "distance_bp", "classification"], rows)


def cmd_utr(args: argparse.Namespace) -> None:
    a, b = read_tsv(args.source_a), load_model_groups(args.source_b)
    rows = []
    for tx in a:
        partners = b.get(tx["original_gene_id"], [])
        partner = partners[0] if partners else None
        available_a, available_b = bool(tx.get("utr_coordinates")), bool(partner and partner.get("utr_coordinates"))
        same = available_a and available_b and tx["utr_coordinates"] == partner["utr_coordinates"]
        rows.append([args.pair_id, tx["original_transcript_id"], partner["original_transcript_id"] if partner else "", int(available_a), int(available_b), int(same), "missing" if not (available_a and available_b) else ("concordant" if same else "discordant")])
    write_tsv(args.output, ["pair_id", "source_a_transcript_id", "source_b_transcript_id", "source_a_utr_available", "source_b_utr_available", "exact_match", "classification"], rows)


def cmd_summary(args: argparse.Namespace) -> None:
    chain_rows = read_tsv(args.intron)
    end_rows = read_tsv(args.ends)
    utr_rows = read_tsv(args.utr)
    summary = [[args.pair_id, len(chain_rows), sum(r.get("classification") == "exact" for r in chain_rows), sum(r.get("classification") == "partial" for r in chain_rows), len(end_rows), sum(r.get("classification") == "within_threshold" for r in end_rows), len(utr_rows), sum(r.get("classification") == "concordant" for r in utr_rows)]]
    write_tsv(args.output, ["pair_id", "n_intron_comparisons", "n_exact_intron_chains", "n_partial_intron_chains", "n_end_comparisons", "n_ends_within_threshold", "n_utr_comparisons", "n_exact_utr_matches"], summary)
    disagreements = [[args.pair_id, r.get("source_a_transcript_id", ""), r.get("source_b_transcript_id", ""), "intron_chain", r.get("classification", "")] for r in chain_rows if r.get("classification") == "disagree"]
    disagreements += [[args.pair_id, r.get("source_a_transcript_id", ""), r.get("source_b_transcript_id", ""), "transcript_end", r.get("classification", "")] for r in end_rows if r.get("classification") == "outside_threshold"]
    write_tsv(args.disagreements, ["pair_id", "source_a_transcript_id", "source_b_transcript_id", "feature", "classification"], disagreements)


def cmd_evidence(args: argparse.Namespace) -> None:
    evidence = {row.get("transcript_id"): row for row in read_tsv(args.evidence)}
    rows = []
    for row in read_tsv(args.models):
        ev = evidence.get(row.get("original_transcript_id"), {})
        rows.append([row.get("pair_id", ""), row.get("source", ""), row.get("original_transcript_id", ""), ev.get("read_count", ""), ev.get("junction_support", ""), ev.get("full_length_support", ""), ev.get("source_tool", "")])
    write_tsv(args.output, ["pair_id", "source", "transcript_id", "read_count", "junction_support", "full_length_support", "source_tool"], rows)


def cmd_bed12gtf(args: argparse.Namespace) -> None:
    """Convert BED12 to GTF with explicitly synthetic gene identifiers."""
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
            tx_attrs = f'gene_id "{gene}"; transcript_id "{tx}"; source_format "BED12";'
            out.write(f"{chrom}\tbed12_to_gtf\ttranscript\t{int(start0)+1}\t{end}\t{score}\t{strand or '.'}\t.\t{tx_attrs}\n")
            for size, offset in zip(sizes, offsets):
                start = int(start0) + offset + 1
                block_end = start + size - 1
                out.write(f"{chrom}\tbed12_to_gtf\texon\t{start}\t{block_end}\t{score}\t{strand or '.'}\t.\t{tx_attrs}\n")


def cmd_models2gtf(args: argparse.Namespace) -> None:
    """Create a GffCompare-compatible GTF from canonical models."""
    with open(args.input, encoding="utf-8") as handle, open(args.output, "w", encoding="utf-8") as out:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row.get("transcript_start") or not row.get("transcript_end"):
                continue
            tx = row["original_transcript_id"]
            gene = row["original_gene_id"]
            if row.get("source_format") == "bed12":
                gene = f"synthetic_gene:{gene}"
            tags = f'gene_id "{gene}"; transcript_id "{tx}";'
            out.write(f'{row["contig"]}\tpairwise_qc\ttranscript\t{row["transcript_start"]}\t{row["transcript_end"]}\t.\t{row["strand"]}\t.\t{tags}\n')
            for exon in filter(None, row.get("exon_coordinates", "").split(";")):
                start, end = exon.split("-", 1)
                out.write(f'{row["contig"]}\tpairwise_qc\texon\t{start}\t{end}\t.\t{row["strand"]}\t.\t{tags}\n')


def cmd_status(args: argparse.Namespace) -> None:
    write_tsv(args.output, ["pair_id", "metric", "status", "reason"], [[args.pair_id, metric, status, reason] for metric, status, reason in args.entry])


def cmd_report(args: argparse.Namespace) -> None:
    summary = read_tsv(args.summary)
    row = summary[0] if summary else {}
    Path(args.output).write_text("# Pairwise annotation QC report\n\n" + f"- Pair: `{args.pair_id}`\n- Intron-chain comparisons: {row.get('n_intron_comparisons', 0)}\n- Exact intron chains: {row.get('n_exact_intron_chains', 0)}\n- End comparisons: {row.get('n_end_comparisons', 0)}\n- UTR comparisons: {row.get('n_utr_comparisons', 0)}\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    v = sub.add_parser("validate"); v.add_argument("--pair-id", required=True); v.add_argument("--assembly", required=True); v.add_argument("--source-a", required=True); v.add_argument("--source-b", required=True); v.add_argument("--source-a-path", required=True); v.add_argument("--source-b-path", required=True); v.add_argument("--source-a-format", required=True); v.add_argument("--source-b-format", required=True); v.add_argument("--source-a-assembly", default=""); v.add_argument("--source-b-assembly", default=""); v.add_argument("--output", required=True); v.set_defaults(func=cmd_validate)
    a = sub.add_parser("audit"); a.add_argument("--models", required=True); a.add_argument("--pair-id", required=True); a.add_argument("--source", required=True); a.add_argument("--source-format", required=True); a.add_argument("--conversion-loss", default="0"); a.add_argument("--output", required=True); a.add_argument("--id-map", required=True); a.set_defaults(func=cmd_audit)
    n = sub.add_parser("normalise"); n.add_argument("--input", required=True); n.add_argument("--source", required=True); n.add_argument("--pair-id", required=True); n.add_argument("--assembly", required=True); n.add_argument("--format", required=True); n.add_argument("--output", required=True); n.set_defaults(func=lambda x: cmd_normalise(x))
    for name, fn in (("intron", cmd_intron), ("ends", cmd_ends), ("utr", cmd_utr)):
        q = sub.add_parser(name); q.add_argument("--source-a", required=True); q.add_argument("--source-b", required=True); q.add_argument("--pair-id", required=True); q.add_argument("--output", required=True); q.set_defaults(func=fn)
    sub.choices["ends"].add_argument("--threshold", type=int, default=50)
    s = sub.add_parser("summary"); s.add_argument("--pair-id", required=True); s.add_argument("--intron", required=True); s.add_argument("--ends", required=True); s.add_argument("--utr", required=True); s.add_argument("--output", required=True); s.add_argument("--disagreements", required=True); s.set_defaults(func=cmd_summary)
    e = sub.add_parser("evidence"); e.add_argument("--models", required=True); e.add_argument("--evidence", required=True); e.add_argument("--output", required=True); e.set_defaults(func=cmd_evidence)
    b = sub.add_parser("bed12gtf"); b.add_argument("--input", required=True); b.add_argument("--output", required=True); b.set_defaults(func=cmd_bed12gtf)
    m = sub.add_parser("models2gtf"); m.add_argument("--input", required=True); m.add_argument("--output", required=True); m.set_defaults(func=cmd_models2gtf)
    st = sub.add_parser("status"); st.add_argument("--pair-id", required=True); st.add_argument("--output", required=True); st.add_argument("--entry", nargs=3, action="append", required=True); st.set_defaults(func=cmd_status)
    r = sub.add_parser("report"); r.add_argument("--pair-id", required=True); r.add_argument("--summary", required=True); r.add_argument("--output", required=True); r.set_defaults(func=cmd_report)
    return p


def cmd_normalise(args: argparse.Namespace) -> None:
    rows, warnings = parse_models(args.input, args.source, args.pair_id, args.assembly, args.format)
    write_tsv(args.output, MODEL_HEADER, rows)
    if warnings:
        Path(args.output + ".warnings").write_text("\n".join(warnings) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parsed = parser().parse_args()
    parsed.func(parsed)
