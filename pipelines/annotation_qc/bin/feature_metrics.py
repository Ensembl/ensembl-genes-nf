#!/usr/bin/env python3
"""Calculate feature-level annotation metrics from A's normalized TSV."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from statistics import mean, median

from qc_utils import parse_biotypes, read_normalized, resolve_selected_features, write_metrics


PHASE_LABELS = {"0": "3n", "1": "3n+1", "2": "3n+2"}
EXON_BINS = ((1, 9), (10, 24), (25, 49), (50, 99), (100, 199), (200, 499), (500, 999), (1000, None))
INTRON_BINS = ((1, 19), (20, 99), (100, 499), (500, 999), (1000, 4999), (5000, 9999), (10000, 49999), (50000, None))


def summary(values):
    if not values:
        return None
    ordered = sorted(values)
    return len(ordered), mean(ordered), median(ordered), ordered[0], ordered[-1]


def add_summary(rows, name, values, suffix="", include_count=False):
    result = summary(values)
    if result is None:
        if include_count:
            rows.append((f"{name} count", 0))
        rows.extend([(f"{name} mean{suffix}", "N/A"), (f"{name} median{suffix}", "N/A"), (f"{name} min{suffix}", "N/A"), (f"{name} max{suffix}", "N/A")])
        return
    count, avg, med, minimum, maximum = result
    if include_count:
        rows.append((f"{name} count", count))
    rows.extend([(f"{name} mean{suffix}", f"{avg:.2f}"), (f"{name} median{suffix}", f"{med:.2f}"), (f"{name} min{suffix}", f"{minimum:.2f}"), (f"{name} max{suffix}", f"{maximum:.2f}")])


def count_pct(count, total):
    return f"{count} (N/A)" if not total else f"{count} ({count / total * 100:.2f}%)"


def histogram(values, bins):
    result = []
    for minimum, maximum in bins:
        if maximum is None:
            label = f">={minimum} bp"
            count = sum(value >= minimum for value in values)
        else:
            label = f"{minimum}-{maximum} bp"
            count = sum(minimum <= value <= maximum for value in values)
        result.append((label, count))
    return result


def merged_exons(exons):
    merged = []
    for start, end in sorted((feature.start, feature.end) for feature in exons):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def transcript_intron_lengths(exons_by_transcript):
    lengths = []
    for exons in exons_by_transcript.values():
        by_locus = defaultdict(list)
        for exon in exons:
            by_locus[(exon.seqid, exon.strand)].append(exon)
        for locus_exons in by_locus.values():
            merged = merged_exons(locus_exons)
            lengths.extend(merged[index][0] - merged[index - 1][1] for index in range(1, len(merged)))
    return [length for length in lengths if length > 0]


def utr_type(feature_type):
    value = feature_type.lower().replace("-", "_").replace("'", "")
    if value in {"five_prime_utr", "five_prime_untranslated_region", "5_prime_utr", "5_utr", "5utr"}:
        return "five_prime"
    if value in {"three_prime_utr", "three_prime_untranslated_region", "3_prime_utr", "3_utr", "3utr"}:
        return "three_prime"
    return None


def utr_summary(segments_by_transcript, coding_transcripts):
    lengths, counts, present = [], [], set()
    for transcript_id in coding_transcripts:
        segments = segments_by_transcript.get(transcript_id, [])
        if segments:
            present.add(transcript_id)
            lengths.append(sum(item.end - item.start for item in segments))
            counts.append(len(segments))
    return present, lengths, counts


def collect_flagged(exons_by_transcript, cds_without_phase, cds_lengths, min_intron, min_exon):
    rows = []
    for transcript_id, length in sorted(cds_lengths.items()):
        if length % 3:
            rows.append(("frame_inconsistent_transcript", transcript_id, "transcript", ".", ".", ".", ".", length, f"CDS length modulo 3 = {length % 3}"))
    for transcript_id, exons in sorted(exons_by_transcript.items()):
        by_locus = defaultdict(list)
        for exon in exons:
            length = exon.end - exon.start
            if length < min_exon:
                rows.append(("short_exon", transcript_id, "exon", exon.seqid, exon.start + 1, exon.end, exon.strand, length, f"exon length < {min_exon} bp"))
            by_locus[(exon.seqid, exon.strand)].append(exon)
        for (seqid, strand), locus_exons in by_locus.items():
            merged = merged_exons(locus_exons)
            for index in range(1, len(merged)):
                start, end = merged[index - 1][1], merged[index][0]
                length = end - start
                if 0 < length < min_intron:
                    rows.append(("short_intron", transcript_id, "intron", seqid, start + 1, end, strand, length, f"intron length < {min_intron} bp"))
    rows.extend(cds_without_phase)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("normalized")
    parser.add_argument("--biotype", action="append")
    parser.add_argument("--exclude-biotype", action="append")
    parser.add_argument("--canonical-only", action="store_true")
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--flagged-features")
    parser.add_argument("--min-intron-length", type=int, default=20)
    parser.add_argument("--min-exon-length", type=int, default=3)
    args = parser.parse_args()

    features = resolve_selected_features(read_normalized(args.normalized), parse_biotypes(args.biotype), parse_biotypes(args.exclude_biotype), args.canonical_only)
    exons = [feature for feature in features if feature.feature_type.lower() == "exon"]
    cds = [feature for feature in features if feature.feature_type.lower() == "cds"]
    exons_by_transcript = defaultdict(list)
    cds_lengths = Counter()
    utrs = {"five_prime": defaultdict(list), "three_prime": defaultdict(list)}
    exon_phases, cds_phases = Counter(), Counter()
    exon_missing = cds_missing = 0
    cds_without_phase = []

    for feature in exons:
        transcript_ids = feature.transcript_ids or feature.parent_ids
        for transcript_id in transcript_ids:
            exons_by_transcript[transcript_id].append(feature)
        if feature.phase in {"0", "1", "2"}:
            exon_phases[feature.phase] += 1
        else:
            exon_missing += 1

    for feature in cds:
        transcript_ids = feature.transcript_ids or feature.parent_ids
        if feature.phase in {"0", "1", "2"}:
            cds_phases[feature.phase] += 1
        else:
            cds_missing += 1
            cds_without_phase.append(("cds_without_usable_phase", ",".join(transcript_ids) or ".", "CDS", feature.seqid, feature.start + 1, feature.end, feature.strand, feature.end - feature.start, "phase is not 0, 1, or 2"))
        for transcript_id in transcript_ids:
            cds_lengths[transcript_id] += feature.end - feature.start

    for feature in features:
        kind = utr_type(feature.feature_type)
        if kind:
            for transcript_id in feature.transcript_ids or feature.parent_ids:
                utrs[kind][transcript_id].append(feature)

    exon_lengths = [feature.end - feature.start for feature in exons]
    intron_lengths = transcript_intron_lengths(exons_by_transcript)
    phase_counts = exon_phases if exon_phases else cds_phases
    phase_missing = exon_missing if exon_phases else cds_missing
    phase_prefix = "Exon" if exon_phases else "CDS"
    phase_source = "exon phase column" if exon_phases else "CDS phase column"
    phase_total = sum(phase_counts.values())
    frame_consistent = sum(length % 3 == 0 for length in cds_lengths.values())
    frame_inconsistent = len(cds_lengths) - frame_consistent
    coding_transcripts = set(cds_lengths)
    rows = [("Exon count", len(exons)), ("Transcripts with exon records", len(exons_by_transcript)), ("Coding transcripts with CDS records", len(cds_lengths)), ("Phase source", phase_source if phase_total else "N/A"), (f"{phase_prefix} records with usable phase", phase_total)]
    for phase, label in PHASE_LABELS.items():
        rows.append((f"{phase_prefix} phase {label}", count_pct(phase_counts[phase], phase_total)))
    rows.append((f"{phase_prefix} records without usable phase", phase_missing if phase_total else "N/A"))
    add_summary(rows, "Exon length", exon_lengths, " (bp)")
    rows.extend((f"Exon length {label}", count_pct(count, len(exon_lengths))) for label, count in histogram(exon_lengths, EXON_BINS))
    add_summary(rows, "Intron length", intron_lengths, " (bp)")
    rows.extend((f"Intron length {label}", count_pct(count, len(intron_lengths))) for label, count in histogram(intron_lengths, INTRON_BINS))

    five, three = set(), set()
    for kind, label in (("five_prime", "5' UTR"), ("three_prime", "3' UTR")):
        present, lengths, counts = utr_summary(utrs[kind], coding_transcripts)
        if kind == "five_prime":
            five = present
        else:
            three = present
        rows.append((f"Coding transcripts with {label}", count_pct(len(present), len(coding_transcripts))))
        add_summary(rows, f"{label} length per transcript", lengths, " (bp)")
        add_summary(rows, f"{label} exon count per transcript", counts)
    rows.extend([( "Coding transcripts with both UTRs", count_pct(len(five & three), len(coding_transcripts))), ("Coding transcripts with neither UTR", count_pct(len(coding_transcripts - five - three), len(coding_transcripts))), ("Frame-consistent coding transcripts", count_pct(frame_consistent, len(cds_lengths))), ("Frame-inconsistent coding transcripts", count_pct(frame_inconsistent, len(cds_lengths)))])

    if args.flagged_features:
        flagged = collect_flagged(exons_by_transcript, cds_without_phase, cds_lengths, args.min_intron_length, args.min_exon_length)
        with open(args.flagged_features, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("flag_type", "transcript_id", "feature_type", "seqid", "start", "end", "strand", "length", "details"))
            writer.writerows(flagged)
        rows.extend([("Flagged records written", len(flagged)), ("Short intron threshold (bp)", f"< {args.min_intron_length}"), ("Short exon threshold (bp)", f"< {args.min_exon_length}")])
    write_metrics(rows, args.output_csv)


if __name__ == "__main__":
    main()

