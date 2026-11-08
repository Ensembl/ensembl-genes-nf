#!/usr/bin/env python3
"""Attach a small, auditable evolutionary-history record to each instance.

The input evidence is deliberately flat and read-only: one JSONL/TSV row per
species candidate. This stage does not select an ortholog, calculate a score,
or change the branch-length denominator.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


STATUS_FIELDS = (
    "orthology_status", "synteny_status", "cds_status",
    "leader_status", "orf_status", "translation_status", "path_status",
)

STATUS_ALIASES = {
    "uorf_status": "orf_status",
    "orf_state": "orf_status",
    "translation_state": "translation_status",
}


def read_rows(path):
    with open(path) as handle:
        first = next((line for line in handle if line.strip()), "")
    with open(path) as handle:
        if first.startswith("{"):
            yield from (json.loads(line) for line in handle if line.strip())
        else:
            yield from csv.DictReader(handle, delimiter="\t")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def candidate(row):
    species = row.get("species")
    if not species:
        raise ValueError("candidate evidence row is missing species")
    fields = dict(row)
    for source, target in STATUS_ALIASES.items():
        if not fields.get(target) and fields.get(source):
            fields[target] = fields[source]
    result = {
        "species": species,
        "reference_gene": fields.get("reference_gene") or None,
        "reference_transcript": fields.get("reference_transcript") or None,
        "reference_cds": fields.get("reference_cds") or None,
        "reference_uorf": fields.get("reference_uorf") or None,
        "candidate_gene": fields.get("candidate_gene") or fields.get("gene_id") or None,
        "candidate_transcript": fields.get("candidate_transcript") or fields.get("transcript_id") or None,
        "locus": fields.get("locus") or None,
        "orthology_source": fields.get("orthology_source") or None,
        "provenance": fields.get("provenance") or None,
    }
    for field in STATUS_FIELDS:
        result[field] = fields.get(field) or "unknown"
    if result["path_status"] == "unknown":
        coordinate_fields = [fields.get(key) for key in ("chrom", "start", "end", "strand")]
        if all(value not in (None, "") for value in coordinate_fields):
            try:
                valid_coordinates = int(fields.get("start")) < int(fields.get("end"))
            except (TypeError, ValueError):
                valid_coordinates = False
            valid_strand = fields.get("strand") in ("+", "-", ".")
            result["path_status"] = "valid" if valid_coordinates and valid_strand else "invalid"
    path = {key: fields.get(key) for key in ("source", "chrom", "start", "end", "strand")
            if fields.get(key) not in (None, "")}
    if path:
        result["path"] = path
    if fields.get("confidence"):
        result["confidence"] = fields["confidence"]
    return result


def history_for(instance, candidates, evidence_path):
    by_species = defaultdict(list)
    for item in candidates:
        by_species[item["species"]].append(item)
    species = [
        {"species": name, "candidates": values, "candidate_count": len(values)}
        for name, values in sorted(by_species.items())
    ]
    summary = {
        "species_with_candidates": len(species),
        "candidate_count": len(candidates),
        "species_with_one_to_one": sum(
            any(x["orthology_status"] == "one_to_one" for x in values)
            for values in by_species.values()
        ),
        "species_with_ambiguous_orthology": sum(
            any(x["orthology_status"] in {"one_to_many", "uncertain"} for x in values)
            for values in by_species.values()
        ),
        "species_states": {
            name: species_state(values) for name, values in sorted(by_species.items())
        },
    }
    if candidates:
        reason = "candidate evidence enumerated; no evolutionary score calculated"
        candidate_state = "enumerated"
    else:
        reason = "no candidate evidence supplied for this instance"
        candidate_state = "not_supplied"
    return {
        "state": "uninformative",
        "reason": reason,
        "candidate_state": candidate_state,
        "reference": {
            "instance_id": instance.get("instance_id"),
            "interval_id": instance.get("interval_id"),
            "gene_id": instance.get("gene_id") or instance.get("gene"),
            "transcript_id": instance.get("transcript_id"),
            "cds": instance.get("cds") or instance.get("cds_id"),
            "uorf": instance.get("uorf") or (instance.get("interval_id") if instance.get("class") == "uORF" else None),
            "class": instance.get("class"),
            "chrom": instance.get("chrom"),
            "start": instance.get("start"),
            "end": instance.get("end"),
            "strand": instance.get("strand"),
            "frame": instance.get("frame"),
        },
        "species": species,
        "summary": summary,
        "denominator": {
            "policy": "full_species_tree",
            "changed": False,
            "note": "candidate lookup does not alter ORBL or BLS denominators",
        },
        "provenance": {
            "candidate_evidence": Path(evidence_path).name,
            "candidate_evidence_sha256": sha256(evidence_path),
            "method": "read-only-candidate-evidence-v1",
        },
    }


def species_state(candidates):
    """Summarise candidate states without selecting a path or inventing absence."""
    genes = {x["candidate_gene"] for x in candidates if x.get("candidate_gene")}
    transcripts = {x["candidate_transcript"] for x in candidates if x.get("candidate_transcript")}
    orthology = {x["orthology_status"] for x in candidates}
    synteny = {x["synteny_status"] for x in candidates}
    cds = {x["cds_status"] for x in candidates}
    leaders = {x["leader_status"] for x in candidates}
    orfs = {x["orf_status"] for x in candidates}
    translations = {x["translation_status"] for x in candidates}
    paths = {x["path_status"] for x in candidates}

    def identity_state(values, statuses):
        if "none" in statuses and not values:
            return "absent"
        if len(values) > 1 or statuses & {"one_to_many", "uncertain"}:
            return "ambiguous"
        return "present" if values else "unknown"

    def evidence_state(values, positive, negative):
        has_positive = bool(values & positive)
        has_negative = bool(values & negative)
        has_ambiguous = bool(values & {"ambiguous", "uncertain", "partial", "discordant"})
        if (has_positive and has_negative) or has_ambiguous:
            return "ambiguous"
        if has_positive:
            return "present"
        if has_negative and values <= negative:
            return "absent"
        return "unknown"

    return {
        "gene": identity_state(genes, orthology),
        "transcript": "ambiguous" if len(transcripts) > 1 or orthology & {"one_to_many", "uncertain"}
        else ("present" if transcripts else "unknown"),
        "synteny": evidence_state(synteny, {"supported", "aligned", "present"}, {"discordant", "missing", "none"}),
        "path": evidence_state(paths, {"valid", "continuous", "resolved"}, {"invalid", "missing", "none"}),
        "cds": evidence_state(cds, {"aligned", "complete", "present"}, {"missing", "none"}),
        "leader": evidence_state(leaders, {"aligned", "complete", "present"}, {"missing", "none"}),
        "orf": evidence_state(orfs, {"candidate", "observed", "present"}, {"not_observed", "missing", "none"}),
        "translation": evidence_state(translations, {"observed", "translated", "present"}, {"not_observed", "missing", "none"}),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    evidence_by_instance = defaultdict(list)
    evidence_by_interval = defaultdict(list)
    for row in read_rows(args.evidence):
        item = candidate(row)
        if row.get("instance_id"):
            evidence_by_instance[row["instance_id"]].append(item)
        elif row.get("interval_id"):
            evidence_by_interval[row["interval_id"]].append(item)
        else:
            raise ValueError("candidate evidence row needs instance_id or interval_id")

    with open(args.output, "w") as handle:
        for instance in read_rows(args.instances):
            candidates = list(evidence_by_instance.get(instance.get("instance_id"), []))
            candidates.extend(evidence_by_interval.get(instance.get("interval_id"), []))
            record = dict(instance)
            record["axes"] = dict(instance.get("axes", {}))
            record["axes"]["evolutionary_history"] = history_for(instance, candidates, args.evidence)
            handle.write(json.dumps(record, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
