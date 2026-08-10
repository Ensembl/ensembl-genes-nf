#!/usr/bin/env python3
"""Write transcript models in genePred and BED12 formats from a GTF."""
import argparse
import re
from collections import defaultdict

ATTR = re.compile(r'([A-Za-z0-9_.-]+)\s+"([^"]+)"')


def parse_attributes(value):
    return dict(ATTR.findall(value))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genepred", required=True)
    parser.add_argument("--bed12", required=True)
    args = parser.parse_args()
    transcripts = defaultdict(lambda: {"chrom": None, "strand": "+", "exons": [], "cds": []})
    with open(args.gtf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() not in {"exon", "cds"}:
                continue
            attrs = parse_attributes(fields[8])
            tid = attrs.get("transcript_id") or attrs.get("transcript")
            if not tid:
                continue
            record = transcripts[tid]
            record["chrom"], record["strand"] = fields[0], fields[6]
            start, end = int(fields[3]) - 1, int(fields[4])
            if fields[2].lower() == "exon":
                record["exons"].append((start, end))
            else:
                record["cds"].append((start, end))
    with open(args.genepred, "w") as gp, open(args.bed12, "w") as bed:
        for tid in sorted(transcripts):
            rec = transcripts[tid]
            exons = sorted(rec["exons"])
            if not exons:
                continue
            chrom, strand = rec["chrom"], rec["strand"]
            tx_start, tx_end = exons[0][0], exons[-1][1]
            cds_start = min((x[0] for x in rec["cds"]), default=tx_start)
            cds_end = max((x[1] for x in rec["cds"]), default=tx_start)
            starts = ",".join(str(x[0]) for x in exons) + ","
            ends = ",".join(str(x[1]) for x in exons) + ","
            gp.write("\t".join(map(str, [tid, chrom, strand, tx_start, tx_end,
                cds_start, cds_end, len(exons), starts, ends])) + "\n")
            sizes = ",".join(str(x[1] - x[0]) for x in exons) + ","
            offsets = ",".join(str(x[0] - tx_start) for x in exons) + ","
            bed.write("\t".join(map(str, [chrom, tx_start, tx_end, tid, 0, strand,
                cds_start, cds_end, "0", len(exons), sizes, offsets])) + "\n")


if __name__ == "__main__":
    main()
