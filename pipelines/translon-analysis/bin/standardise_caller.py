#!/usr/bin/env python3
"""Convert caller-native ORF output into the translon-consensus contract."""
import argparse
import csv
import json
import math
import re
from pathlib import Path

FIELDS = ["sample_id", "tool", "chrom", "start", "end", "strand", "frame",
          "transcript_id", "orf_id", "score", "pval", "qval", "extra_json"]

PRICE_LOCATION_RE = re.compile(r"^(?P<chrom>.+?)(?P<strand>[+-]):(?P<blocks>.+)$")


def score(value, pvalue=""):
    try:
        number = float(pvalue)
        if 0 <= number <= 1:
            return max(0, min(1000, round(1000 * (1 - number))))
    except (TypeError, ValueError):
        pass
    try:
        return max(0, min(1000, round(1000 * (1 - math.exp(-abs(float(value)))))))
    except (TypeError, ValueError, OverflowError):
        return 0


def rows_from_raw(raw, tool, sample):
    raw = Path(raw)
    def candidates(*patterns):
        if raw.is_file():
            return [raw] if any(raw.match(pattern) for pattern in patterns) else []
        return sorted({path for pattern in patterns for path in raw.glob(pattern)})

    if tool == "orfquant":
        orfquant_paths = candidates("*_Detected_ORFs.gtf")
        for path in orfquant_paths:
            grouped = {}
            with path.open() as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9 or fields[2].lower() != "cds":
                        continue
                    attrs = {}
                    for item in fields[8].split(";"):
                        item = item.strip()
                        if not item:
                            continue
                        bits = item.split(None, 1)
                        if len(bits) == 2:
                            attrs[bits[0]] = bits[1].strip().strip('"')
                    orf_id = attrs.get("ORF_id") or attrs.get("orf_id") or attrs.get("ID")
                    if not orf_id:
                        continue
                    record = grouped.setdefault(orf_id, {"fields": fields, "attrs": attrs, "blocks": []})
                    record["blocks"].append((int(fields[3]) - 1, int(fields[4])))
            for orf_id, record in grouped.items():
                fields, attrs, blocks = record["fields"], record["attrs"], sorted(record["blocks"])
                start, end = min(x[0] for x in blocks), max(x[1] for x in blocks)
                extra = dict(attrs)
                extra["_blocks"] = blocks
                yield [sample, tool, fields[0], start, end, fields[6], fields[7],
                       attrs.get("transcript_id", ""), orf_id,
                       attrs.get("ORF_pct_P_sites") or attrs.get("ORFs_pM") or fields[5],
                       attrs.get("pval", ""), attrs.get("qval", ""), json.dumps(extra)]
        if orfquant_paths:
            return
    if tool in {"rpbp", "iribo", "price", "riborf"}:
        if tool == "price":
            for path in candidates("*.tsv"):
                with path.open() as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        orf_id = (row.get("Id") or "").strip()
                        match = PRICE_LOCATION_RE.match((row.get("Location") or "").strip())
                        if not orf_id or not match:
                            continue
                        blocks = []
                        for token in match.group("blocks").split("|"):
                            try:
                                left, right = token.split("-", 1)
                                blocks.append((int(left), int(right)))
                            except (TypeError, ValueError):
                                continue
                        if not blocks:
                            continue
                        blocks.sort()
                        orf_type = row.get("Type") or row.get("ORF_type") or ""
                        pval = row.get("p value") or row.get("p_value") or ""
                        transcript_id = orf_id
                        if orf_type:
                            match_id = re.match(rf"^(.+)_{re.escape(orf_type)}_\d+$", orf_id)
                            if match_id:
                                transcript_id = match_id.group(1)
                        extra = dict(row)
                        extra["_blocks"] = blocks
                        yield [sample, tool, match.group("chrom"), blocks[0][0], blocks[-1][1],
                               match.group("strand"), ".", transcript_id, orf_id,
                               pval, pval, "", json.dumps(extra)]
            return
        for bed in candidates("*.bed"):
            with bed.open() as handle:
                for line in handle:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    chrom, start, end = fields[:3]
                    yield [sample, tool, chrom, start, end,
                           fields[4] if len(fields) > 4 and fields[4] in {"+", "-"} else ".",
                           None, None, fields[3] if len(fields) > 3 else "orf",
                           fields[5] if len(fields) > 5 else "", None, None, "{}"]
        return

    if tool == "ribotricer":
        # Ribotricer encodes transcript-relative ORF coordinates in ORF_ID
        # and reports translating/non-translating candidates in a dedicated
        # table rather than a generic start/end table.
        candidates = sorted(raw.glob("*translating_ORFs.tsv"))
        for path in candidates:
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    if row.get("status") != "translating":
                        continue
                    parts = row.get("ORF_ID", "").rsplit("_", 3)
                    if len(parts) != 4:
                        continue
                    transcript_id, start, end, _length = parts
                    try:
                        # Ribotricer ORF_ID coordinates are one-based and
                        # inclusive; the consensus contract is BED-like.
                        start, end = int(start) - 1, int(end)
                    except ValueError:
                        continue
                    extra = dict(row)
                    yield [sample, tool, row.get("chrom") or ".", start, end,
                           row.get("strand") or ".", ".",
                           row.get("transcript_id") or transcript_id,
                           row.get("ORF_ID") or "orf",
                           row.get("phase_score") or row.get("read_density") or "",
                           "", "", json.dumps(extra)]
        return

    if tool == "ribotish":
        candidates = sorted(raw.glob("*.txt"))
        for path in candidates:
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    if not (row.get("RiboPStatus") or "").startswith("T"):
                        continue
                    genome_pos = row.get("GenomePos", "")
                    try:
                        chrom, coords, strand = genome_pos.rsplit(":", 2)
                        genomic_start, genomic_end = [int(x) for x in coords.split("-", 1)]
                    except ValueError:
                        continue
                    blocks = []
                    for block in (row.get("Blocks") or "").split(","):
                        try:
                            block_start, block_end = [int(x) for x in block.split("-", 1)]
                        except ValueError:
                            continue
                        blocks.append([block_start - 1, block_end])
                    if not blocks:
                        blocks = [[genomic_start - 1, genomic_end]]
                    extra = dict(row)
                    extra["_blocks"] = blocks
                    yield [sample, tool, chrom, genomic_start - 1, genomic_end,
                           strand, ".", row.get("Tid") or "", row.get("Tid") or "orf",
                           row.get("RiboPvalue") or "", row.get("TISPvalue") or "",
                           row.get("FisherQvalue") or "", json.dumps(extra)]
        return

    for path in candidates("*.tsv", "*.txt"):
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                start, end = row.get("start"), row.get("end")
                if not start or not end:
                    continue
                yield [sample, tool, row.get("chrom") or row.get("chr") or ".", start, end,
                       row.get("strand") or ".", row.get("frame") or ".",
                       row.get("transcript_id") or row.get("transcript") or "",
                       row.get("orf_id") or row.get("orf") or row.get("name") or "orf",
                       row.get("score") or row.get("TPM") or row.get("reads") or "",
                       row.get("pval") or "", row.get("qval") or "",
                       json.dumps({k: v for k, v in row.items()
                                   if k not in {"start", "end", "chrom", "chr", "strand", "frame",
                                                "transcript_id", "transcript", "orf_id", "orf", "name",
                                                "score", "TPM", "reads", "pval", "qval"}})]


def load_transcripts(gtf):
    transcripts = {}
    if not gtf:
        return transcripts
    with Path(gtf).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2].lower() not in {"exon", "cds"}:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            parent = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("transcript")
            if not parent:
                continue
            parent = parent.replace("transcript:", "").split(",")[0]
            record = transcripts.setdefault(parent, {"chrom": fields[0], "strand": fields[6], "exons": []})
            start, end = int(fields[3]) - 1, int(fields[4])
            if fields[2].lower() == "exon":
                record["exons"].append((start, end))
    for record in transcripts.values():
        record["exons"].sort(reverse=record["strand"] == "-")
    return transcripts


def transcript_blocks(row, transcript):
    start, end = int(row[3]), int(row[4])
    cursor = 0
    blocks = []
    for exon_start, exon_end in transcript["exons"]:
        exon_len = exon_end - exon_start
        left, right = max(start, cursor), min(end, cursor + exon_len)
        if left < right:
            if transcript["strand"] == "+":
                blocks.append((exon_start + left - cursor, exon_start + right - cursor))
            else:
                blocks.append((exon_end - (right - cursor), exon_end - (left - cursor)))
        cursor += exon_len
    return sorted(blocks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bed12", required=True, type=Path)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--adapter-version", default="1")
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows_from_raw(args.raw, args.tool, args.sample))
    transcripts = load_transcripts(args.gtf)
    blocks_by_row = {}
    for index, row in enumerate(rows):
        transcript_id = row[7] or ""
        try:
            native_blocks = json.loads(row[12]).get("_blocks")
            if native_blocks:
                blocks_by_row[index] = [tuple(block) for block in native_blocks]
        except (TypeError, ValueError, json.JSONDecodeError):
            pass
        transcript = transcripts.get(transcript_id.split(".")[0])
        if row[2] in {"", ".", "None"} and transcript:
            blocks = transcript_blocks(row, transcript)
            if blocks:
                row[2], row[3], row[4], row[5] = transcript["chrom"], min(x[0] for x in blocks), max(x[1] for x in blocks), transcript["strand"]
                blocks_by_row[index] = blocks
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(FIELDS)
        writer.writerows(rows)

    with args.bed12.open("w") as handle:
        for index, row in enumerate(rows):
            chrom, start, end = row[2], int(row[3]), int(row[4])
            if end <= start:
                continue
            strand = row[5] if row[5] in {"+", "-"} else "."
            name = f"{row[8] or 'orf'}|{row[7] or 'tx'}|{args.tool}"
            blocks = blocks_by_row.get(index, [(start, end)])
            sizes = ",".join(str(block_end - block_start) for block_start, block_end in blocks) + ","
            offsets = ",".join(str(block_start - min(x[0] for x in blocks)) for block_start, _ in blocks) + ","
            chrom_start, chrom_end = min(x[0] for x in blocks), max(x[1] for x in blocks)
            handle.write("\t".join(map(str, [chrom, chrom_start, chrom_end, name, score(row[9], row[10]),
                strand, chrom_start, chrom_end, "0,0,0", len(blocks), sizes, offsets])) + "\n")


if __name__ == "__main__":
    main()
