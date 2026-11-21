#!/usr/bin/env python3
"""Export a versioned per-sample processing bundle from a local RDP release."""

import argparse
import csv
import hashlib
import json
import re
import sqlite3
from collections import Counter, defaultdict, deque
from contextlib import closing, contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory


SETUP = Path(__file__).resolve().parents[1] / "species_run_setup"
RUN = re.compile(r"[SED]RR[0-9]+")
PROJECT = re.compile(r"PRJ(?:NA|EB|DB)[0-9]+")
SLUG = re.compile(r"[a-z0-9_]+")
ENSEMBL_FTP = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"
COMMUNITY_PROVIDERS = {
    "GCA_001433935.1": "community_rapdb",
    "GCA_902167145.1": "community_cshl",
}
SAMPLE_COLUMNS = """
Run BioProject_id SRAStudy Experiment BioSample TaxID ScientificName LIBRARYTYPE
LibraryStrategy LibraryLayout LibrarySource INHIBITOR Nuclease FRACTION TISSUE
CELL_LINE CONDITION Strain UMI Adapter spots bases avgLength INSDC_status
process_status verified
""".split()
REFERENCE_COLUMNS = """
species slug taxids config_status beta_assembly beta_provider beta_geneset_date
genome_fasta_url gtf_url reference_choice_note
""".split()
PROCESSING_PARAMS = {
    "run_matrix_mode": False,
    "enable_unique_reads_tracking": False,
    "run_organism_setup": False,
    "auto_build_indices": False,
    "fetch": True,
    "force_fetch": True,
    "exact_ena_fastq": True,
    "rpf_extraction_method": "getrpf",
    "getrpf_command": "extract",
    "getrpf_method": "structure",
    "run_rrna_filter": False,
}


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest()}


def check_snapshot(db):
    for suffix in ("-wal", "-journal"):
        if Path(str(db) + suffix).exists():
            raise ValueError(f"Use a closed SQLite release snapshot. Found {db.name}{suffix}")


def migrate_ensembl_urls(ref):
    accession = ref["beta_assembly"]
    legacy = re.compile(rf"{re.escape(ENSEMBL_FTP)}/[^/]+/{re.escape(accession)}/")
    match = re.fullmatch(r"(GC[AF])_([0-9]{3})([0-9]{3})([0-9]{3})\.([0-9]+)", accession)
    if not match:
        raise ValueError(f"Invalid assembly accession {accession!r}")
    if ref["beta_provider"] == "community":
        ref["beta_provider"] = COMMUNITY_PROVIDERS.get(accession, "community")
    release = "/".join((*match.groups(), ref["beta_provider"], ref["beta_geneset_date"]))
    root = f"{ENSEMBL_FTP}/{release}"
    for field, directory in (("genome_fasta_url", "genome"), ("gtf_url", "geneset")):
        url = ref[field]
        if legacy.match(url):
            filename = url.rsplit("/", 1)[-1]
            if directory == "genome" and filename.endswith(".fa.gz"):
                filename = filename[:-3] + ".bgz"
            ref[field] = f"{root}/{directory}/{filename}"


def read_references(path):
    references = {}
    by_taxon = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ref = {key: row[key].strip() for key in REFERENCE_COLUMNS}
            slug = ref["slug"]
            if not SLUG.fullmatch(slug) or slug in references:
                raise ValueError(f"Invalid or repeated reference slug {slug!r}")
            ref["reference_id"] = ""
            if ref["config_status"] in {"beta_url_config", "pre_release_url_config"}:
                if not all(ref[key] for key in REFERENCE_COLUMNS[4:9]):
                    raise ValueError(f"Incomplete Beta reference for {slug}")
                migrate_ensembl_urls(ref)
                prefix = "pre_release" if ref["beta_provider"] == "pre_release" else "beta"
                version = "{prefix}_{beta_geneset_date}_{beta_assembly}_{beta_provider}".format(
                    prefix=prefix, **ref
                )
                ref["reference_id"] = f"{slug}__{version}"
            taxids = [item.partition(":")[0] for item in ref["taxids"].split(";")]
            if any(not taxid.isdecimal() or int(taxid) < 1 for taxid in taxids):
                raise ValueError(f"Invalid TaxID list for {slug}")
            ref["taxids"] = ";".join(taxids)
            references[slug] = ref
            for taxid in taxids:
                taxon = (taxid, ref["species"].casefold())
                if taxon in by_taxon:
                    raise ValueError(f"Repeated reference mapping for {taxon}")
                by_taxon[taxon] = ref
    return references, by_taxon


def read_review_queue(path):
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not {"Run", "reasons"}.issubset(reader.fieldnames or []):
            raise ValueError("The RDP review queue must contain Run and reasons columns")
        return {row["Run"].strip(): row["reasons"].strip() or "RDP review queue" for row in reader}


def read_inventory(db, by_taxon, review_queue):
    inventory = []
    seen = set()
    columns = ", ".join(SAMPLE_COLUMNS)
    with closing(sqlite3.connect(db.as_uri() + "?mode=ro", uri=True)) as conn:
        conn.row_factory = sqlite3.Row
        for row in conn.execute(f"SELECT {columns} FROM main_sample ORDER BY Run"):
            run = {key: "" if value is None else str(value).strip() for key, value in dict(row).items()}
            accession = run["Run"]
            if not accession or accession in seen:
                raise ValueError(f"Missing or repeated run accession {accession!r} in {db}")
            seen.add(accession)

            ref = by_taxon.get((run["TaxID"], run["ScientificName"].casefold()), {})
            run["reference_id"] = ref.get("reference_id", "")
            run["species_slug"] = ref.get("slug", "")
            run["reference_status"] = (
                "configured" if run["reference_id"] else ref.get("config_status", "unmapped_taxon")
            )
            run["study_accession"] = run["BioProject_id"]
            run["rdp_review_reason"] = review_queue.get(accession, "")
            reasons = []
            if not RUN.fullmatch(accession):
                reasons.append("invalid_run_accession")
            if not PROJECT.fullmatch(run["study_accession"]):
                reasons.append("unresolved_bioproject")
            if not run["LIBRARYTYPE"]:
                reasons.append("unknown_library_type")
            elif run["LIBRARYTYPE"] != "Ribo-Seq":
                reasons.append("outside_baseline_library_type")
            if run["LibraryLayout"] != "SINGLE":
                reasons.append("outside_single_end_baseline")
            if run["INSDC_status"].casefold() not in ("", "public"):
                reasons.append("archive_status_needs_review")
            if run["rdp_review_reason"]:
                reasons.append("rdp_metadata_review")
            if not run["reference_id"]:
                reasons.append("reference_pending")
            run["processing_status"] = "held" if reasons else "candidate"
            run["hold_reasons"] = ";".join(reasons)
            # RDP metadata alone does not resolve run-level FASTQs or checksums.
            run["raw_input_status"] = "unresolved"
            inventory.append(run)
    if not inventory:
        raise ValueError(f"No sample rows found in {db}")
    if unknown := review_queue.keys() - seen:
        raise ValueError(f"Review queue contains runs outside this snapshot: {', '.join(sorted(unknown)[:5])}")
    return inventory


def stratum(run):
    inhibitor = run["INHIBITOR"].casefold() or "unknown"
    nuclease = run["Nuclease"].casefold() or "unknown"
    length = int(run["avgLength"] or 0)
    spots = int(run["spots"] or 0)
    length_bin = "unknown" if length <= 0 else "lt35" if length < 35 else "35-74" if length < 75 else "75plus"
    depth_bin = "unknown" if spots <= 0 else "lt10m" if spots < 10_000_000 else "10m-49m" if spots < 50_000_000 else "50mplus"
    return (inhibitor, nuclease, length_bin, depth_bin)


def study_balanced(runs, seed):
    def rank(value):
        return hashlib.sha256(f"{seed}:{value}".encode()).digest()

    studies = defaultdict(lambda: defaultdict(deque))
    for run in sorted(runs, key=lambda row: rank(row["Run"])):
        studies[run["study_accession"]][stratum(run)].append(run)
    queues = {
        study: deque(sorted(groups, key=lambda key: rank((study, key))))
        for study, groups in studies.items()
    }
    while queues:
        groups = defaultdict(deque)
        for study, strata in queues.items():
            groups[strata[0]].append(study)
        # Each study supplies at most one run per round. Spread that round
        # across the observed inhibitor, nuclease, depth and length strata.
        for key in groups:
            groups[key] = deque(sorted(groups[key], key=rank))
        active = deque(sorted(groups, key=rank))
        while active:
            key = active.popleft()
            study = groups[key].popleft()
            if groups[key]:
                active.append(key)
            queue = queues[study]
            queue.popleft()
            yield studies[study][key].popleft()
            if studies[study][key]:
                queue.append(key)
            if not queue:
                del queues[study]


def choose_pilot(plan, inventory, candidates):
    by_run = {row["Run"]: row for row in inventory}
    controls = plan["controls"]
    quotas = plan["quotas"]
    invalid_control = [acc for acc in controls if acc not in by_run
                       or by_run[acc]["processing_status"] != "candidate"
                       or by_run[acc]["species_slug"] not in quotas]
    if invalid_control:
        raise ValueError(f"Pilot controls missing, held or outside pilot taxa: {', '.join(invalid_control)}")

    selection = []
    for slug, count in quotas.items():
        pool = candidates.get(slug, [])
        named = [row for row in pool if row["Run"] in controls]
        if type(count) is not int or count <= 0 or not len(named) <= count <= len(pool):
            raise ValueError(f"Pilot quota for {slug} cannot be filled from {len(pool)} candidates")
        general = study_balanced((row for row in pool if row["Run"] not in controls), plan["seed"])
        selected = named + [next(general) for _ in range(count - len(named))]
        for row in sorted(selected, key=lambda run: run["Run"]):
            intent = controls.get(row["Run"], "")
            selection.append({
                **{key: row[key] for key in ("Run", "study_accession", "species_slug", "reference_id")},
                "selection_group": "named_control" if intent else "study_balanced",
                "control_intent": intent,
                **dict(zip(("inhibitor", "nuclease", "length_bin", "depth_bin"), stratum(row))),
            })
    return selection


def write_table(path, fields, data, *, delimiter="\t", lineterminator="\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter=delimiter, extrasaction="ignore",
                                lineterminator=lineterminator)
        writer.writeheader()
        writer.writerows(data)


def write_sample_sheet(directory, slug, rows):
    write_table(directory / f"{slug}_samplesheet.csv", ["Run", "study_accession"], rows,
                delimiter=",", lineterminator="\r\n")


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_bundle(out, inventory, references, candidates, pilot, sources):
    write_table(out / "run_inventory.tsv", list(inventory[0]), inventory)
    write_table(out / "references.tsv", ["reference_id", *REFERENCE_COLUMNS], references.values())
    for slug, rows in sorted(candidates.items()):
        write_sample_sheet(out / "sample_sheets", slug, rows)
    for slug, ref in references.items():
        if not ref["reference_id"]:
            continue
        write_json(out / "organism_setup" / f"{slug}.json", {
            "organism": ref["species"],
            "ensembl_version": ref["reference_id"].split("__", 1)[1],
            "genome_fasta_url": ref["genome_fasta_url"],
            "genome_gtf_url": ref["gtf_url"],
            "run_organism_setup": True,
            "auto_build_indices": False,
            "download_method": "url",
            "rrna_source": "gtf",
        })
    write_json(out / "processing.json", PROCESSING_PARAMS)
    summary = defaultdict(Counter)
    for row in inventory:
        key = tuple(row[field] for field in ("TaxID", "ScientificName", "species_slug", "reference_status"))
        summary[key][row["processing_status"]] += 1
    summary_fields = ["TaxID", "ScientificName", "species_slug", "reference_status"]
    write_table(out / "species_summary.tsv", summary_fields + ["total", "candidate", "held"], (
        {**dict(zip(summary_fields, key)), "total": sum(counts.values()),
         "candidate": counts["candidate"], "held": counts["held"]}
        for key, counts in sorted(summary.items())
    ))
    if pilot:
        write_table(out / "pilot" / "selection.tsv", list(pilot[0]), pilot)
        for slug in sorted({row["species_slug"] for row in pilot}):
            rows = [row for row in pilot if row["species_slug"] == slug]
            write_sample_sheet(out / "pilot", slug, rows)
            controls = [row for row in rows if row["selection_group"] == "named_control"]
            if controls:
                write_sample_sheet(out / "controls", slug, controls)
    write_json(out / "provenance.json", {
        "bundle_version": 2,
        "sources": sources,
        "counts": {
            "runs": len(inventory),
            "candidates": sum(map(len, candidates.values())),
            "configured_references": sum(bool(ref["reference_id"]) for ref in references.values()),
            "pilot_runs": len(pilot),
        },
        "policy": {
            "library_type": "Ribo-Seq",
            "layout": "SINGLE",
            "reference_match": "exact TaxID and case-insensitive scientific name in reference manifest",
            "pilot_status": "candidate selection awaiting control validation" if pilot else None,
        },
    })


@contextmanager
def new_bundle(output):
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise ValueError(f"Output already exists. Use a new versioned directory: {output}")
    # Publish a complete bundle only after all inputs and quotas validate.
    with TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as staging:
        bundle = Path(staging) / "bundle"
        bundle.mkdir()
        yield bundle
        bundle.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True, help="Closed RDP release SQLite snapshot")
    parser.add_argument("--references", type=Path, default=SETUP / "species_manifest.tsv")
    parser.add_argument("--review-queue", type=Path, help="Optional review_queue.csv from the same RDP release")
    parser.add_argument("--pilot-plan", type=Path, help="Optional JSON with seed, species quotas and named controls")
    parser.add_argument("--out-dir", type=Path, required=True, help="New directory for a versioned export")
    args = parser.parse_args()
    try:
        db = args.db.resolve(strict=True)
        check_snapshot(db)
        inputs = {
            "database": db,
            "references": args.references,
            "exporter": Path(__file__),
            "review_queue": args.review_queue,
            "pilot_plan": args.pilot_plan,
        }
        sources = {key: fingerprint(path.resolve(strict=True)) for key, path in inputs.items() if path is not None}
        references, by_taxon = read_references(args.references)
        inventory = read_inventory(db, by_taxon, read_review_queue(args.review_queue))
        check_snapshot(db)
        candidates = defaultdict(list)
        for row in inventory:
            if row["processing_status"] == "candidate":
                candidates[row["species_slug"]].append(row)

        pilot = []
        if args.pilot_plan:
            plan = json.loads(args.pilot_plan.read_text(encoding="utf-8"))
            pilot = choose_pilot(plan, inventory, candidates)
            sources["pilot_plan"]["selection_inputs"] = plan
        for key, path in inputs.items():
            if path is not None and fingerprint(path.resolve(strict=True))["sha256"] != sources[key]["sha256"]:
                raise ValueError(f"Input changed during export: {path}")
        with new_bundle(args.out_dir.absolute()) as out:
            write_bundle(out, inventory, references, candidates, pilot, sources)
    except (OSError, KeyError, TypeError, ValueError, sqlite3.Error) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")
    print(f"Exported {len(inventory):,} RDP runs, {sum(map(len, candidates.values())):,} candidates "
          f"and {len(pilot):,} pilot candidates to {args.out_dir}")


if __name__ == "__main__":
    main()
