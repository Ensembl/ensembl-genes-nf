#!/usr/bin/env python3
"""Join validated BED12 models to run-scoped direct alignment observations."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import mmap
import re
import sqlite3
from collections import Counter
from pathlib import Path

from bam_to_features import gzip_text, table

SCHEMA = "model-event-support-v1"


def digest(path: Path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


class FastaReference:
    def __init__(self, path: Path):
        self.handle = path.open("rb")
        self.data = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)
        self.index = {}
        name = None
        offset = width = stride = length = 0
        previous_short = False
        with path.open("rb") as stream:
            while True:
                line_start = stream.tell()
                line = stream.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    if name is not None:
                        self.index[name] = (offset, width, stride, length)
                    name = line[1:].split()[0].decode()
                    offset = width = stride = length = 0
                    previous_short = False
                    continue
                bases = line.rstrip(b"\r\n")
                if not bases:
                    continue
                if not width:
                    offset, width, stride = line_start, len(bases), len(line)
                elif previous_short or len(bases) > width or (len(bases) == width and len(line) != stride):
                    raise ValueError(f"irregular FASTA wrapping: {path}")
                previous_short = len(bases) < width
                length += len(bases)
        if name is not None:
            self.index[name] = (offset, width, stride, length)
        if not self.index:
            raise ValueError(f"empty FASTA: {path}")

    def fetch(self, contig: str, start: int, end: int):
        if contig not in self.index:
            return None
        offset, width, stride, length = self.index[contig]
        if start < 0 or end > length:
            return None
        return bytes(self.data[offset + (position // width) * stride + position % width] for position in range(start, end)).decode().upper()

    def motif(self, contig: str, start: int, end: int, strand: str):
        first, last = self.fetch(contig, start, start + 2), self.fetch(contig, end - 2, end)
        if not first or not last or strand not in ("+", "-"):
            return "not_available"
        if strand == "-":
            complement = str.maketrans("ACGTN", "TGCAN")
            first, last = last.translate(complement)[::-1], first.translate(complement)[::-1]
        motif = first + "-" + last
        return motif + (":canonical" if motif in {"GT-AG", "GC-AG", "AT-AC"} else ":noncanonical")

    def close(self):
        self.data.close()
        self.handle.close()


def rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parts(chain: str):
    return chain.split(";") if chain else []


def contiguous_subset(shorter, longer):
    size = len(shorter)
    return bool(size) and any(longer[i:i + size] == shorter for i in range(len(longer) - size + 1))


def parse_bed(path: Path):
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                raise ValueError(f"invalid validated BED12: {path}")
            left, right = int(fields[1]), int(fields[2])
            sizes = [int(v) for v in fields[10].rstrip(",").split(",")]
            starts = [left + int(v) for v in fields[11].rstrip(",").split(",")]
            if len(sizes) != int(fields[9]) or len(starts) != len(sizes):
                raise ValueError(f"invalid BED12 blocks: {path}")
            introns = [(starts[i] + sizes[i], starts[i + 1]) for i in range(len(starts) - 1)]
            yield fields[0], left, right, fields[3], fields[5], introns


def create_database(path: Path, support_root: Path):
    db = sqlite3.connect(path)
    db.execute("CREATE TABLE runs (run_id TEXT PRIMARY KEY, assembly_id TEXT, sample_id TEXT, library_id TEXT, technology TEXT, min_mapq INTEGER, read_set_id TEXT)")
    db.execute("CREATE TABLE chains (run_id TEXT, contig TEXT, chain TEXT, records INTEGER)")
    db.execute("CREATE TABLE chain_member (run_id TEXT, contig TEXT, intron TEXT, chain TEXT, records INTEGER)")
    db.execute("CREATE TABLE junctions (run_id TEXT, contig TEXT, start INTEGER, end INTEGER, records INTEGER)")
    db.execute("CREATE TABLE ends (run_id TEXT, contig TEXT, chain TEXT, left_end INTEGER, right_end INTEGER, records INTEGER)")
    support_dirs = sorted(p for p in support_root.iterdir() if p.is_dir())
    if not support_dirs:
        raise ValueError("no direct-support run directories were staged")
    for directory in support_dirs:
        metadata = json.loads((directory / "evidence_metadata.json").read_text())
        if metadata.get("schema") != "direct-alignment-support-v1":
            raise ValueError(f"unsupported direct evidence schema: {directory}")
        run_id = metadata["run_id"]
        min_mapq = int(metadata["support_filter"].rsplit("-q ", 1)[1])
        db.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)", (run_id, metadata["assembly_id"], metadata["sample_id"], metadata["library_id"], metadata["technology"], min_mapq, metadata["read_set_id"]))
        for row in rows(directory / "chain_observations.tsv.gz"):
            chain, count = row["intron_chain"], int(row["alignment_records"])
            db.execute("INSERT INTO chains VALUES (?,?,?,?)", (run_id, row["contig"], chain, count))
            db.executemany("INSERT INTO chain_member VALUES (?,?,?,?,?)", ((run_id, row["contig"], intron, chain, count) for intron in parts(chain)))
        for row in rows(directory / "junction_observations.tsv.gz"):
            db.execute("INSERT INTO junctions VALUES (?,?,?,?,?)", (run_id, row["contig"], int(row["intron_start"]), int(row["intron_end"]), int(row["alignment_records"])))
        for row in rows(directory / "end_observations.tsv.gz"):
            db.execute("INSERT INTO ends VALUES (?,?,?,?,?,?)", (run_id, row["contig"], row["intron_chain"], int(row["left_aligned_end"]), int(row["right_aligned_end"]), int(row["alignment_records"])))
    db.execute("CREATE INDEX chain_key ON chains(run_id,contig,chain)")
    db.execute("CREATE INDEX member_key ON chain_member(run_id,contig,intron)")
    db.execute("CREATE INDEX junction_key ON junctions(run_id,contig,start,end)")
    db.execute("CREATE INDEX end_key ON ends(run_id,contig,chain)")
    db.commit()
    return db


def support_for(db, run_id, contig, chain, introns, left, right):
    if introns:
        exact = db.execute("SELECT COALESCE(SUM(records),0) FROM chains WHERE run_id=? AND contig=? AND chain=?", (run_id, contig, chain)).fetchone()[0]
        candidates = {}
        for intron in parts(chain):
            for observed, count in db.execute("SELECT chain,SUM(records) FROM chain_member WHERE run_id=? AND contig=? AND intron=? GROUP BY chain", (run_id, contig, intron)):
                candidates[observed] = count
        partial = sum(count for observed, count in candidates.items() if observed != chain and contiguous_subset(parts(observed), parts(chain)))
        junction_counts = [db.execute("SELECT COALESCE(SUM(records),0) FROM junctions WHERE run_id=? AND contig=? AND start=? AND end=?", (run_id, contig, start, end)).fetchone()[0] for start, end in introns]
        end_where = "run_id=? AND contig=? AND chain=?"
        end_params = (run_id, contig, chain)
        relation = "exact_chain"
    else:
        exact, partial, junction_counts = 0, 0, []
        end_where = "run_id=? AND contig=? AND chain='' AND left_end<=? AND right_end>=?"
        end_params = (run_id, contig, right - 1, left)
        relation = "single_exon_overlap"
    depth = db.execute("SELECT COALESCE(SUM(records),0) FROM ends WHERE " + end_where, end_params).fetchone()[0]
    return exact, partial, min(junction_counts) if junction_counts else "NA", end_where, end_params, depth, relation


def endpoint_counts(db, end_where, end_params, coordinate, side):
    column = "left_end" if side == "left" else "right_end"
    query = "SELECT COALESCE(SUM(records),0)," + ",".join(
        f"COALESCE(SUM(CASE WHEN ABS({column}-?)<={tolerance} THEN records ELSE 0 END),0)"
        for tolerance in (0, 5, 10, 20)
    ) + " FROM ends WHERE " + end_where
    depth, exact, near5, near10, near20 = db.execute(query, (coordinate, coordinate, coordinate, coordinate, *end_params)).fetchone()
    if not depth:
        return depth, exact, near5, near10, near20, "NA", "NA", "NA", "NA"
    mode, mode_count = db.execute(
        f"SELECT {column},SUM(records) FROM ends WHERE {end_where} GROUP BY {column} ORDER BY SUM(records) DESC,{column} LIMIT 1",
        end_params,
    ).fetchone()
    mean_distance = db.execute(
        f"SELECT SUM(records*ABS({column}-?))*1.0/SUM(records) FROM ends WHERE {end_where}",
        (coordinate, *end_params),
    ).fetchone()[0]
    return depth, exact, near5, near10, near20, mode, mode_count, abs(mode - coordinate), mean_distance


def event_rows(db):
    query = "SELECT assembly_id,contig,strand,event_type,coordinate,backend,run_id,sample_id,technology,library_id,model_id,depth,near5,near10,near20,relation,read_set_id,exact FROM event_rows ORDER BY assembly_id,contig,strand,event_type,coordinate,backend,run_id,model_id"
    yield from db.execute(query)


def depth_bin(depth):
    if depth == 0:
        return "0"
    if depth < 5:
        return "1-4"
    if depth < 10:
        return "5-9"
    if depth < 20:
        return "10-19"
    if depth < 50:
        return "20-49"
    return "50+"


def chain_bin(count):
    if count == 0:
        return "0"
    if count == 1:
        return "1"
    if count < 5:
        return "2-4"
    if count < 10:
        return "5-9"
    return "10+"


def write_cluster_sensitivity(output_dir):
    source = output_dir / "event_support.tsv.gz"
    target = output_dir / "event_cluster_sensitivity.tsv.gz"
    fields = ("schema_version", "assembly_id", "contig", "strand", "event_type", "tolerance_bp", "cluster_start", "cluster_end", "exact_event_count", "best_raw_depth", "any_raw_support", "any_recurrent_support", "any_shared_caller_support")
    with gzip_text(target) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for tolerance in (5, 10, 20):
            current = None
            def flush():
                if current is not None:
                    writer.writerow((SCHEMA, *current["key"], tolerance, current["start"], current["end"], current["count"], current["depth"], int(current["raw"]), int(current["recurrent"]), int(current["shared"])))
            for row in rows(source):
                key = (row["assembly_id"], row["contig"], row["strand"], row["event_type"])
                coordinate = int(row["coordinate"])
                if current is None or key != current["key"] or coordinate - current["end"] > tolerance:
                    flush()
                    current = {"key": key, "start": coordinate, "end": coordinate, "count": 0, "depth": 0, "raw": False, "recurrent": False, "shared": False}
                current["end"] = coordinate
                current["count"] += 1
                current["depth"] = max(current["depth"], int(row["best_raw_depth"]))
                current["raw"] |= int(row["raw_run_count"]) > 0
                current["recurrent"] |= int(row["raw_sample_count"]) >= 2 or int(row["raw_technology_count"]) >= 2
                current["shared"] |= int(row["caller_count"]) >= 2
            flush()


def write_chain_events(db, output_dir):
    fields = ("schema_version", "assembly_id", "contig", "strand", "intron_chain", "intron_count", "model_count", "caller_count", "run_count", "sample_count", "library_count", "read_set_count", "technology_count", "raw_caller_count", "raw_run_count", "raw_sample_count", "raw_library_count", "raw_read_set_count", "raw_technology_count", "best_exact_chain_records", "best_partial_chain_records", "best_minimum_junction_records", "review_tier")
    query = "SELECT assembly_id,contig,strand,intron_chain,model_id,backend,run_id,sample_id,library_id,read_set_id,technology,exact_records,partial_records,minimum_junction_records FROM chain_event_rows ORDER BY assembly_id,contig,strand,intron_chain"

    def output_rows():
        for key, group in itertools.groupby(db.execute(query), lambda row: row[:4]):
            models = set()
            callers = set()
            runs = set()
            samples = set()
            libraries = set()
            read_sets = set()
            technologies = set()
            raw_callers = set()
            raw_runs = set()
            raw_samples = set()
            raw_libraries = set()
            raw_read_sets = set()
            raw_technologies = set()
            best_exact = best_partial = best_junction = 0
            for row in group:
                _, _, _, chain, model, caller, run, sample, library, read_set, technology, exact, partial, junction = row
                models.add(model)
                callers.add(caller)
                runs.add(run)
                samples.add(sample)
                libraries.add(library)
                read_sets.add(read_set)
                technologies.add(technology)
                best_exact = max(best_exact, exact)
                best_partial = max(best_partial, partial)
                best_junction = max(best_junction, junction)
                if exact:
                    raw_callers.add(caller)
                    raw_runs.add(run)
                    raw_samples.add(sample)
                    raw_libraries.add(library)
                    raw_read_sets.add(read_set)
                    raw_technologies.add(technology)
            tier = ("partial_only" if best_exact == 0 and best_partial > 0 else
                    "model_only" if best_exact == 0 else
                    "recurrent_raw_supported" if len(raw_samples) >= 2 or len(raw_technologies) >= 2 else
                    "depth_supported" if best_exact >= 20 else "candidate_raw_supported")
            yield (SCHEMA, *key, chain.count(";") + 1, len(models), len(callers), len(runs), len(samples), len(libraries), len(read_sets), len(technologies), len(raw_callers), len(raw_runs), len(raw_samples), len(raw_libraries), len(raw_read_sets), len(raw_technologies), best_exact, best_partial, best_junction, tier)

    table(output_dir / "chain_event_support.tsv.gz", fields, output_rows())


def write_events(db, output_dir):
    caller_fields = ("schema_version", "assembly_id", "contig", "strand", "event_type", "coordinate", "backend", "run_id", "sample_id", "technology", "library_id", "read_set_id", "model_count", "best_raw_depth", "best_exact_end_records", "best_near10", "best_fraction_20plus")
    exact_fields = ("schema_version", "assembly_id", "contig", "strand", "event_type", "coordinate", "model_count", "caller_count", "run_count", "sample_count", "library_count", "read_set_count", "technology_count", "exact_run_count", "exact_sample_count", "exact_library_count", "exact_read_set_count", "exact_technology_count", "raw_run_count", "raw_sample_count", "raw_library_count", "raw_read_set_count", "raw_technology_count", "best_raw_depth", "best_exact_end_records", "best_near10", "best_fraction_20plus", "review_tier")
    sensitivity = Counter()
    tiers = Counter()
    depths = Counter()
    with gzip_text(output_dir / "event_caller_run_support.tsv.gz") as caller_handle, gzip_text(output_dir / "event_support.tsv.gz") as exact_handle:
        caller_writer = csv.writer(caller_handle, delimiter="\t", lineterminator="\n")
        exact_writer = csv.writer(exact_handle, delimiter="\t", lineterminator="\n")
        caller_writer.writerow(caller_fields)
        exact_writer.writerow(exact_fields)
        key_func = lambda row: row[:5]
        for key, group in itertools.groupby(event_rows(db), key_func):
            model_rows = list(group)
            callers = set()
            runs = set()
            samples = set()
            libraries = set()
            read_sets = set()
            technologies = set()
            raw_runs = set()
            raw_samples = set()
            raw_libraries = set()
            raw_read_sets = set()
            raw_technologies = set()
            exact_runs = set()
            exact_samples = set()
            exact_libraries = set()
            exact_read_sets = set()
            exact_technologies = set()
            best_depth = best_near10 = 0
            best_exact = 0
            best_fraction20 = None
            for caller_key, sub in itertools.groupby(model_rows, lambda row: (row[5], row[6])):
                same = list(sub)
                backend, run_id = caller_key
                callers.add(backend)
                runs.add(run_id)
                samples.update(row[7] for row in same)
                libraries.update(row[9] for row in same)
                read_sets.update(row[16] for row in same)
                technologies.update(row[8] for row in same)
                if any(row[13] > 0 for row in same):
                    raw_runs.add(run_id)
                    raw_samples.add(same[0][7])
                    raw_technologies.add(same[0][8])
                    raw_libraries.add(same[0][9])
                    raw_read_sets.add(same[0][16])
                if any(row[17] > 0 for row in same):
                    exact_runs.add(run_id)
                    exact_samples.add(same[0][7])
                    exact_libraries.add(same[0][9])
                    exact_read_sets.add(same[0][16])
                    exact_technologies.add(same[0][8])
                best = max(same, key=lambda row: (row[11], row[13], row[10]))
                eligible = [row for row in same if row[11] >= 20]
                fraction20 = max((row[13] / row[11] for row in eligible), default=None)
                best_depth = max(best_depth, best[11])
                best_near10 = max(best_near10, *(row[13] for row in same))
                best_exact = max(best_exact, *(row[17] for row in same))
                if fraction20 is not None:
                    best_fraction20 = max(best_fraction20 or 0, fraction20)
                caller_row = (*key, backend, run_id, same[0][7], same[0][8], same[0][9], same[0][16], len(same), best[11], max(row[17] for row in same), max(row[13] for row in same), fraction20 if fraction20 is not None else "NA")
                caller_writer.writerow((SCHEMA, *caller_row))
            recurrent = len(raw_samples) >= 2 or len(raw_technologies) >= 2
            if best_near10 == 0:
                tier = "model_only"
            elif best_depth < 20:
                tier = "candidate_raw_supported"
            elif recurrent and len(callers) >= 2 and best_fraction20 is not None and best_fraction20 >= 0.5:
                tier = "strong_shared_reproducible"
            elif recurrent:
                tier = "reproducible_supported"
            else:
                tier = "supported_single_observation"
            exact_writer.writerow((SCHEMA, *key, len({row[10] for row in model_rows}), len(callers), len(runs), len(samples), len(libraries), len(read_sets), len(technologies), len(exact_runs), len(exact_samples), len(exact_libraries), len(exact_read_sets), len(exact_technologies), len(raw_runs), len(raw_samples), len(raw_libraries), len(raw_read_sets), len(raw_technologies), best_depth, best_exact, best_near10, best_fraction20 if best_fraction20 is not None else "NA", tier))
            tiers[(key[3], tier)] += 1
            depths[(key[3], depth_bin(best_depth))] += 1
            for tolerance, column in ((5, 12), (10, 13), (20, 14)):
                tolerance_rows = [row for row in model_rows if row[column] > 0]
                tolerance_recurrent = len({row[7] for row in tolerance_rows}) >= 2 or len({row[8] for row in tolerance_rows}) >= 2
                for threshold in (1, 5, 10, 20, 50):
                    prefix = (key[3], tolerance, threshold)
                    sensitivity[(*prefix, "candidate_events")] += 1
                    if tolerance_rows:
                        sensitivity[(*prefix, "any_raw_support")] += 1
                    if any(row[11] >= threshold and row[column] > 0 for row in model_rows):
                        sensitivity[(*prefix, "depth_threshold_reached")] += 1
                    if tolerance_recurrent:
                        sensitivity[(*prefix, "recurrent_raw_support")] += 1
                    if tolerance_recurrent and any(row[11] >= threshold and row[column] > 0 for row in model_rows):
                        sensitivity[(*prefix, "depth_and_recurrent")] += 1
                    eligible_fractions = [row[column] / row[11] for row in model_rows if row[11] >= threshold]
                    if tolerance_recurrent and len(callers) >= 2 and eligible_fractions and max(eligible_fractions) >= 0.5:
                        sensitivity[(*prefix, "shared_recurrent_fraction_ge_0.5")] += 1
    stages = ("candidate_events", "any_raw_support", "depth_threshold_reached", "recurrent_raw_support", "depth_and_recurrent", "shared_recurrent_fraction_ge_0.5")
    sensitivity_rows = ((SCHEMA, event_type, tolerance, threshold, stage, sensitivity[(event_type, tolerance, threshold, stage)]) for event_type in ("5prime", "3prime") for tolerance in (5, 10, 20) for threshold in (1, 5, 10, 20, 50) for stage in stages)
    table(output_dir / "event_support_sensitivity.tsv.gz", ("schema_version", "event_type", "tolerance_bp", "depth_threshold", "stage", "event_count"), sensitivity_rows)
    table(output_dir / "event_tier_summary.tsv", ("schema_version", "event_type", "review_tier", "exact_event_count"), ((SCHEMA, *key, value) for key, value in sorted(tiers.items())))
    table(output_dir / "event_best_depth_summary.tsv", ("schema_version", "event_type", "best_raw_depth_bin", "exact_event_count"), ((SCHEMA, *key, value) for key, value in sorted(depths.items())))
    write_cluster_sensitivity(output_dir)
    db.commit()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models-dir", required=True, type=Path)
    parser.add_argument("--model-records", required=True, type=Path)
    parser.add_argument("--direct-support-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    db_path = args.output_dir / "join.sqlite"
    db_path.unlink(missing_ok=True)
    db = create_database(db_path, args.direct_support_dir)
    reference = FastaReference(args.reference)
    db.execute("CREATE TABLE event_rows (assembly_id TEXT, contig TEXT, strand TEXT, event_type TEXT, coordinate INTEGER, backend TEXT, run_id TEXT, sample_id TEXT, technology TEXT, library_id TEXT, model_id TEXT, depth INTEGER, near5 INTEGER, near10 INTEGER, near20 INTEGER, relation TEXT, read_set_id TEXT, exact INTEGER)")
    db.execute("CREATE TABLE chain_event_rows (assembly_id TEXT, contig TEXT, strand TEXT, intron_chain TEXT, model_id TEXT, backend TEXT, run_id TEXT, sample_id TEXT, library_id TEXT, read_set_id TEXT, technology TEXT, exact_records INTEGER, partial_records INTEGER, minimum_junction_records INTEGER)")
    runs = list(db.execute("SELECT run_id,assembly_id,sample_id,library_id,technology,min_mapq,read_set_id FROM runs ORDER BY run_id"))
    records = json.loads(args.model_records.read_text())
    model_paths = {path.name: path for path in args.models_dir.iterdir() if path.is_file()}
    index_fields = ("schema_version", "model_id", "normalized_id", "native_transcript_id", "model_artifact_sha256", "backend", "scope", "dataset_id", "assembly_id", "contig", "strand", "start", "end", "intron_chain", "intron_count", "model_origin")
    support_fields = ("schema_version", "model_id", "run_id", "read_set_id", "sample_id", "library_id", "technology", "assembly_id", "backend", "scope", "contig", "strand", "intron_chain", "exact_chain_records", "partial_chain_records", "minimum_junction_records", "terminal_depth", "support_status", "terminal_relation", "construction_relation", "full_length_status", "molecule_count", "filter_min_mapq")
    event_fields = ("schema_version", "model_id", "run_id", "read_set_id", "sample_id", "library_id", "technology", "backend", "scope", "assembly_id", "contig", "strand", "event_type", "coordinate", "terminal_relation", "raw_chain_depth", "end_records_exact", "end_records_within_5bp", "end_records_within_10bp", "end_records_within_20bp", "fraction_within_10bp", "observed_mode_coordinate", "observed_mode_records", "distance_to_mode_bp", "mean_absolute_end_distance_bp", "support_status", "construction_relation", "polya_support_status")
    junction_fields = ("schema_version", "model_id", "run_id", "read_set_id", "assembly_id", "contig", "strand", "intron_start", "intron_end", "chain_position", "alignment_records", "splice_motif", "support_status")
    chain_depth = Counter()
    with gzip_text(args.output_dir / "model_index.tsv.gz") as index_handle, gzip_text(args.output_dir / "model_support.tsv.gz") as support_handle, gzip_text(args.output_dir / "model_event_support.tsv.gz") as event_handle, gzip_text(args.output_dir / "model_junction_support.tsv.gz") as junction_handle:
        index_writer = csv.writer(index_handle, delimiter="\t", lineterminator="\n")
        support_writer = csv.writer(support_handle, delimiter="\t", lineterminator="\n")
        event_writer = csv.writer(event_handle, delimiter="\t", lineterminator="\n")
        junction_writer = csv.writer(junction_handle, delimiter="\t", lineterminator="\n")
        index_writer.writerow(index_fields)
        support_writer.writerow(support_fields)
        event_writer.writerow(event_fields)
        junction_writer.writerow(junction_fields)
        for record in sorted(records, key=lambda value: (value["backend"], value["scope"], value["bed12_name"])):
            path = model_paths[record["bed12_name"]]
            artifact_sha256 = digest(path)
            backend, scope = record["backend"], record["scope"]
            dataset_id = record.get("run_id") if scope == "run" else record.get("cohort_id")
            selected_runs = [run for run in runs if scope == "cohort" or run[0] == record.get("run_id")]
            if not selected_runs:
                raise ValueError(f"model has no matching run evidence: {record}")
            if len({run[1] for run in selected_runs}) != 1:
                raise ValueError(f"model spans multiple assemblies: {record}")
            for contig, left, right, normalized_id, strand, introns in parse_bed(path):
                chain = ";".join(f"{start}-{end}" for start, end in introns)
                model_id = hashlib.sha256(f"{artifact_sha256}|{backend}|{scope}|{dataset_id}|{normalized_id}".encode()).hexdigest()
                prefix = f"{backend}.{scope}."
                match = re.match(r"^" + re.escape(prefix) + r"model_\d+\.(.*)$", normalized_id)
                native_id = match.group(1) if match else "unknown"
                index_writer.writerow((SCHEMA, model_id, normalized_id, native_id, artifact_sha256, backend, scope, dataset_id, selected_runs[0][1], contig, strand, left, right, chain, len(introns), "caller_reconstructed"))
                for run_id, assembly, sample, library, technology, min_mapq, read_set_id in selected_runs:
                    exact, partial, min_junction, end_where, end_params, terminal_depth, relation = support_for(db, run_id, contig, chain, introns, left, right)
                    if strand not in ("+", "-"):
                        status = "ambiguous_strand"
                    elif terminal_depth or exact or partial:
                        status = "supported"
                    else:
                        status = "not_observed"
                    construction = "used_to_build"
                    support_writer.writerow((SCHEMA, model_id, run_id, read_set_id, sample, library, technology, assembly, backend, scope, contig, strand, chain, exact, partial, min_junction, terminal_depth, status, relation, construction, "unknown", "NA", min_mapq))
                    if introns:
                        db.execute("INSERT INTO chain_event_rows VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (assembly, contig, strand, chain, model_id, backend, run_id, sample, library, read_set_id, technology, exact, partial, min_junction))
                    for position, (intron_start, intron_end) in enumerate(introns, 1):
                        junction_count = db.execute("SELECT COALESCE(SUM(records),0) FROM junctions WHERE run_id=? AND contig=? AND start=? AND end=?", (run_id, contig, intron_start, intron_end)).fetchone()[0]
                        junction_writer.writerow((SCHEMA, model_id, run_id, read_set_id, assembly, contig, strand, intron_start, intron_end, position, junction_count, reference.motif(contig, intron_start, intron_end, strand), "supported" if junction_count else "not_observed"))
                    if strand not in ("+", "-"):
                        continue
                    for event_type, coordinate, side in (("5prime", left if strand == "+" else right - 1, "left" if strand == "+" else "right"), ("3prime", right - 1 if strand == "+" else left, "right" if strand == "+" else "left")):
                        depth, exact_end, near5, near10, near20, mode, mode_count, mode_distance, mean_distance = endpoint_counts(db, end_where, end_params, coordinate, side)
                        fraction = near10 / depth if depth else "NA"
                        event_status = "supported" if near10 else ("not_observed" if depth else "no_matching_chain")
                        event_writer.writerow((SCHEMA, model_id, run_id, read_set_id, sample, library, technology, backend, scope, assembly, contig, strand, event_type, coordinate, relation, depth, exact_end, near5, near10, near20, fraction, mode, mode_count, mode_distance, mean_distance, event_status, construction, "not_available"))
                        chain_depth[(backend, event_type, chain_bin(len(introns)), depth_bin(depth))] += 1
                        db.execute("INSERT INTO event_rows VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (assembly, contig, strand, event_type, coordinate, backend, run_id, sample, technology, library, model_id, depth, near5, near10, near20, relation, read_set_id, exact_end))
    db.execute("CREATE INDEX event_key ON event_rows(assembly_id,contig,strand,event_type,coordinate,backend,run_id)")
    write_events(db, args.output_dir)
    write_chain_events(db, args.output_dir)
    table(args.output_dir / "chain_length_depth_composition.tsv.gz", ("schema_version", "backend", "event_type", "intron_count_bin", "raw_depth_bin", "model_run_event_count"), ((SCHEMA, *key, value) for key, value in sorted(chain_depth.items())))
    metadata = {"schema": SCHEMA, "event_key": ["assembly_id", "contig", "strand", "event_type", "coordinate"], "coordinates": "zero_based_aligned_end_base", "default_terminal_tolerance_bp": 10, "provisional_depth_threshold": 20, "molecule_counts": "not_available"}
    (args.output_dir / "support_schema.json").write_text(json.dumps(metadata, sort_keys=True, indent=2) + "\n")
    db.close()
    reference.close()
    db_path.unlink()


if __name__ == "__main__":
    main()
