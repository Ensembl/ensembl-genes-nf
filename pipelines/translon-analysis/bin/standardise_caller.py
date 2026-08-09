#!/usr/bin/env python3
"""Convert caller-native ORF output into the translon-consensus contract."""
import argparse
import csv
import json
import math
from pathlib import Path

FIELDS = ["sample_id", "tool", "chrom", "start", "end", "strand", "frame",
          "transcript_id", "orf_id", "score", "pval", "qval", "extra_json"]


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
    if tool in {"ribotaper", "rpbp"}:
        candidates = sorted(raw.glob("*.bed"))
        for bed in candidates:
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

    candidates = sorted(raw.glob("*.tsv"))
    for path in candidates:
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
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows_from_raw(args.raw, args.tool, args.sample))
    transcripts = load_transcripts(args.gtf)
    blocks_by_row = {}
    for index, row in enumerate(rows):
        transcript_id = row[7] or ""
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
