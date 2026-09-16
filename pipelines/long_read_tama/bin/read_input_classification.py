#!/usr/bin/env python3
"""Pure, offline-friendly read representation inspection and approval gate.

The module intentionally keeps declared platform and observed representation
separate.  Its functions are also used by the command line tools below, which
makes the decision rules straightforward to unit test without archive access.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import shutil
import urllib.parse
import urllib.request
import io
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO

RUN = re.compile(r"^(?:SRR|ERR|DRR)[0-9]+$")
MD5 = re.compile(r"^[0-9a-fA-F]{32}$")
SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
UUID = re.compile(r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[1-5][0-9a-fA-F]{3}-[89abAB][0-9a-fA-F]{3}-[0-9a-fA-F]{12}$")
SUBREAD = re.compile(r"^([^/\s]+)/([0-9]+)/([0-9]+_[0-9]+)$")
CCS = re.compile(r"^([^/\s]+)/([0-9]+)/ccs$")
SAFE_URI = re.compile(r"^(?:https?|ftp)://[^\s]+$", re.I)
SAFE_FILE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

CLASSIFICATION_FIELDS = [
    "run_accession", "tissue", "description", "manifest_declared_platform",
    "ena_platform", "ena_instrument_model", "sra_platform", "sra_spot_group",
    "library_strategy", "library_source", "library_selection",
    "ncbi_original_alias", "ncbi_original_representation",
    "selected_artifact_uri", "selected_artifact_md5", "selected_artifact_basename",
    "submitted_representation", "header_representation", "header_records_sampled",
    "header_subread_count", "header_ccs_count", "header_ont_uuid_count",
    "header_malformed_count", "classification", "confidence", "status", "reason_codes",
    "proposed_action", "review_decision", "reviewer", "reviewed_at",
]
ARTIFACT_FIELDS = ["run_accession", "source", "uri", "basename", "md5", "format", "artifact_role", "selected"]
APPROVED_FIELDS = ["run_accession", "tissue", "description", "classification", "proposed_action",
                   "selected_artifact_uri", "selected_artifact_md5", "selected_artifact_basename",
                   "minimap2_preset", "expected_header_representation", "classification_report_sha256",
                   "reviewer", "reviewed_at"]


def first_token(header: str) -> str:
    return header[1:].split()[0] if header.startswith("@") and header[1:].split() else ""


def header_signature(token: str, header: str = "") -> tuple[str, str | None]:
    match = SUBREAD.fullmatch(token)
    if match:
        return "PACBIO_SUBREAD", f"{match.group(1)}/{match.group(2)}"
    if CCS.fullmatch(token):
        match = CCS.fullmatch(token)
        return "PACBIO_CCS", f"{match.group(1)}/{match.group(2)}"
    if UUID.fullmatch(token) and any(k in header.lower() for k in ("runid=", "ch=", "barcode=")):
        return "ONT", token
    return "UNKNOWN", None


def _records(handle: TextIO) -> Iterable[tuple[str, str, str]]:
    """Yield complete FASTQ records, supporting multiline sequence/quality."""
    while True:
        line = handle.readline()
        if not line:
            return
        if not line.startswith("@"):
            raise ValueError("record does not start with @")
        header = line.rstrip("\r\n")
        sequence = []
        while True:
            line = handle.readline()
            if not line:
                raise ValueError("truncated sequence")
            line = line.rstrip("\r\n")
            if line.startswith("+"):
                break
            sequence.append(line)
        seq = "".join(sequence)
        quality = []
        length = 0
        while length < len(seq):
            line = handle.readline()
            if not line:
                raise ValueError("truncated quality")
            line = line.rstrip("\r\n")
            quality.append(line)
            length += len(line)
        qual = "".join(quality)
        if len(qual) != len(seq):
            raise ValueError("sequence and quality lengths differ")
        yield header, seq, qual


def _probe_remote_bytes(uri: str, max_bytes: int = 2_000_000, timeout: int = 20) -> tuple[bytes, str, str]:
    """Fetch only the beginning of an ENA FASTQ, preferring ranged HTTPS."""
    candidates = [uri]
    if uri.startswith("ftp://ftp.sra.ebi.ac.uk/"):
        candidates.insert(0, "https://ftp.sra.ebi.ac.uk/" + uri.removeprefix("ftp://ftp.sra.ebi.ac.uk/"))
    last_error = ""
    for candidate in candidates:
        try:
            request = urllib.request.Request(candidate, headers={"Range": f"bytes=0-{max_bytes - 1}"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read(max_bytes), candidate, "https-range" if candidate.startswith("https://") else "bounded-stream"
        except (OSError, urllib.error.URLError) as exc:
            last_error = str(exc)
    raise OSError(f"remote FASTQ probe failed for {uri}: {last_error}")


def probe_fastq(path_or_uri: str, limit: int = 100, probe_bytes: int = 2_000_000) -> dict:
    path_or_uri = str(path_or_uri)
    opener = open
    probe_method = "local"
    bytes_read = 0
    probe_uri = path_or_uri
    if re.match(r"^https?://", path_or_uri, re.I) or path_or_uri.startswith("ftp://"):
        data, probe_uri, probe_method = _probe_remote_bytes(path_or_uri, probe_bytes)
        bytes_read = len(data)
        import gzip
        raw = gzip.GzipFile(fileobj=io.BytesIO(data)) if path_or_uri.lower().endswith(".gz") else io.BytesIO(data)
        handle = io.TextIOWrapper(raw)
    else:
        import gzip
        handle = gzip.open(path_or_uri, "rt") if path_or_uri.lower().endswith(".gz") else opener(path_or_uri)
    counts = Counter()
    molecules: set[str] = set()
    headers = []
    sampled = 0
    malformed = 0
    failure = ""
    try:
        iterator = _records(handle)
        while sampled < limit:
            try:
                header, _seq, _qual = next(iterator)
            except StopIteration:
                break
            except (OSError, ValueError) as exc:
                malformed += 1
                failure = str(exc)
                break
            token = first_token(header)
            representation, molecule = header_signature(token, header)
            counts[representation] += 1
            if molecule:
                molecules.add(molecule)
            headers.append(token)
            sampled += 1
    except (OSError, ValueError) as exc:
        malformed += 1
        failure = str(exc)
    finally:
        handle.close()
    decisive = {key for key in ("PACBIO_SUBREAD", "PACBIO_CCS", "ONT") if counts[key]}
    representation = next(iter(decisive)) if len(decisive) == 1 else ("MIXED" if decisive else "UNKNOWN")
    return {"header_representation": representation, "records_sampled": sampled,
            "subread_count": counts["PACBIO_SUBREAD"], "ccs_count": counts["PACBIO_CCS"],
            "ont_uuid_count": counts["ONT"], "malformed_count": malformed,
            "first_tokens": headers, "molecule_count": len(molecules), "duplicate_count": 0, "probe_failure": failure,
            "probe_method": probe_method, "probe_uri": probe_uri, "probe_bytes": bytes_read}


def validate_fastq(path_or_uri: str, expected: str = "UNKNOWN") -> dict:
    """Validate the complete acquired FASTQ and enforce one ID per molecule."""
    path_or_uri = str(path_or_uri)
    import gzip
    handle = gzip.open(path_or_uri, "rt") if path_or_uri.lower().endswith(".gz") else open(path_or_uri)
    read_ids: set[str] = set()
    molecule_ids: set[str] = set()
    duplicates = 0
    records = 0
    representations = set()
    try:
        for header, _seq, _qual in _records(handle):
            token = first_token(header)
            representation, molecule = header_signature(token, header)
            if token in read_ids: duplicates += 1
            read_ids.add(token)
            if molecule:
                if molecule in molecule_ids: duplicates += 1
                molecule_ids.add(molecule)
            representations.add(representation)
            records += 1
    finally:
        handle.close()
    observed = next(iter(representations)) if len(representations) == 1 else ("MIXED" if representations else "UNKNOWN")
    if records == 0: raise ValueError("FASTQ contains no complete records")
    # NCBI can rename a submitted CCS FASTQ and replace its molecule headers
    # while retaining the consensus sequences. Original-format metadata makes
    # UNKNOWN headers acceptable for this explicitly approved representation.
    if expected not in ("UNKNOWN", "NOT_APPLICABLE", "PACBIO_CCS_ORIGINAL", "PACBIO_PROCESSED") and observed != expected: raise ValueError(f"expected {expected}, observed {observed}")
    if duplicates: raise ValueError(f"duplicate read or molecule IDs: {duplicates}")
    return {"records": records, "representation": observed, "distinct_ids": len(read_ids), "distinct_molecules": len(molecule_ids)}


def expand_artifacts(run: str, source: str, uris: str, md5s: str, formats: str = "") -> list[dict]:
    uri_list = [x for x in uris.split(";") if x]
    md5_list = [x for x in md5s.split(";") if x]
    format_list = formats.split(";") if formats else [""] * len(uri_list)
    if len(uri_list) != len(md5_list) or len(format_list) != len(uri_list):
        raise ValueError("semicolon-separated artifact lists have different lengths")
    output = []
    for uri, checksum, fmt in zip(uri_list, md5_list, format_list):
        if not SAFE_URI.fullmatch(uri) or not MD5.fullmatch(checksum):
            raise ValueError(f"invalid artifact URI or checksum for {run}")
        basename = Path(urllib.parse.unquote(urllib.parse.urlparse(uri).path)).name
        lower = basename.lower()
        if not SAFE_FILE.fullmatch(basename):
            raise ValueError(f"unsafe artifact basename: {basename}")
        role = artifact_role(basename)
        output.append({"run_accession": run, "source": source, "uri": uri, "basename": basename,
                       "md5": checksum.lower(), "format": fmt or role, "artifact_role": role, "selected": "false"})
    return output


def artifact_role(basename: str) -> str:
    lower = basename.lower()
    return ("SUBREADS_BAM" if lower.endswith(".subreads.bam") else
            "CCS_BAM" if lower.endswith((".ccs.bam", ".hifi_reads.bam")) else
            "PBI" if lower.endswith(".pbi") else "BAI" if lower.endswith(".bai") else
            "FASTQ" if lower.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")) else "UNKNOWN")


def classify(*, declared_platform: str, ena_platform: str, sra_platform: str,
             library_strategy: str = "", library_source: str = "", library_selection: str = "",
             artifacts: list[dict], probe: dict) -> dict:
    submitted = "UNKNOWN"
    has_subread_bam = any(a["artifact_role"] == "SUBREADS_BAM" for a in artifacts)
    has_ccs_bam = any(a["artifact_role"] == "CCS_BAM" for a in artifacts)
    fastq_artifacts = [a for a in artifacts if a["artifact_role"] == "FASTQ"]
    fastq_subread_hint = any("subread" in a.get("basename", "").lower() for a in fastq_artifacts)
    original_artifacts = [a for a in artifacts if a.get("source") == "NCBI_ORIGINAL"]
    original_ccs_fastq = any(a.get("artifact_role") == "FASTQ" and ".ccs." in a.get("basename", "").lower() for a in original_artifacts)
    original_processed_fastq = any(a.get("artifact_role") == "FASTQ" and a.get("format") == "PACBIO_PROCESSED_FASTQ" for a in original_artifacts)
    original_ccs_bam = any(a.get("artifact_role") == "CCS_BAM" for a in original_artifacts)
    original_subread_bam = any(a.get("artifact_role") == "SUBREADS_BAM" for a in original_artifacts)
    if has_subread_bam or original_subread_bam: submitted = "PACBIO_SUBREAD"
    elif has_ccs_bam or original_ccs_bam or original_ccs_fastq: submitted = "PACBIO_CCS"
    elif original_processed_fastq: submitted = "PACBIO_PROCESSED"
    elif has_ccs_bam: submitted = "PACBIO_CCS"
    elif any(a["artifact_role"] == "FASTQ" for a in artifacts): submitted = probe.get("header_representation", "UNKNOWN")
    header = probe.get("header_representation", "UNKNOWN")
    pacbio_declared = any("PACBIO" in x.upper() or "SMRT" in x.upper() for x in (declared_platform, ena_platform, sra_platform) if x)
    ont_declared = any("ONT" in x.upper() or "OXFORD" in x.upper() for x in (ena_platform, sra_platform) if x)
    reasons = []
    if header == "ONT" and not ont_declared:
        header = "UNKNOWN"
        reasons.append("ONT_UUID_WITHOUT_METADATA_SUPPORT")
    if probe.get("probe_failure"): reasons.append("FASTQ_PROBE_FAILED")
    if probe.get("subread_count"): reasons.append("HEADER_SUBREAD")
    if probe.get("ccs_count"): reasons.append("HEADER_CCS")
    if probe.get("ont_uuid_count"): reasons.append("HEADER_ONT_UUID")
    if pacbio_declared and probe.get("ont_uuid_count"): reasons.append("MANIFEST_HEADER_MISMATCH")
    if ena_platform and sra_platform and ena_platform.lower() != sra_platform.lower() and sra_platform.lower() != "unavailable": reasons.append("ENA_SRA_PLATFORM_MISMATCH")
    if library_strategy: reasons.append(f"LIBRARY_STRATEGY_{library_strategy.upper().replace('-', '_').replace(' ', '_')}")
    if library_source: reasons.append(f"LIBRARY_SOURCE_{library_source.upper().replace('-', '_').replace(' ', '_')}")
    if library_selection: reasons.append(f"LIBRARY_SELECTION_{library_selection.upper().replace('-', '_').replace(' ', '_')}")
    if any(a["artifact_role"] == "SUBREADS_BAM" for a in artifacts): reasons.append("SUBMITTED_SUBREADS_BAM")
    if any(a["artifact_role"] == "CCS_BAM" for a in artifacts): reasons.append("SUBMITTED_CCS_BAM")
    if original_ccs_fastq: reasons.append("NCBI_ORIGINAL_CCS_FASTQ")
    if original_ccs_bam: reasons.append("NCBI_ORIGINAL_CCS_BAM")
    if original_subread_bam: reasons.append("NCBI_ORIGINAL_SUBREADS_BAM")
    if original_processed_fastq: reasons.append("NCBI_ORIGINAL_PROCESSED_FASTQ")
    conflict = (header == "MIXED" or (pacbio_declared and header == "ONT") or (ont_declared and header in ("PACBIO_SUBREAD", "PACBIO_CCS")) or (has_subread_bam and header == "ONT") or (has_ccs_bam and header == "PACBIO_SUBREAD") or ((original_ccs_fastq or original_ccs_bam or original_processed_fastq) and header == "PACBIO_SUBREAD"))
    if conflict: classification, action, confidence = "MIXED" if header == "MIXED" else "CONFLICT", "QUARANTINE", "INSUFFICIENT"
    elif (has_subread_bam or has_ccs_bam) and not pacbio_declared:
        reasons.append("NO_PACBIO_METADATA_SUPPORT")
        classification, action, confidence = "UNKNOWN", "QUARANTINE", "INSUFFICIENT"
    elif has_subread_bam or original_subread_bam: classification, action, confidence = "PACBIO_SUBREAD_BAM", "RUN_CCS_THEN_ALIGN", "CONFIRMED"
    elif has_ccs_bam or original_ccs_bam: classification, action, confidence = "PACBIO_CCS_BAM", "CANONICALISE_CCS_THEN_ALIGN", "CONFIRMED"
    elif original_ccs_fastq: classification, action, confidence = "PACBIO_CCS_FASTQ", "ALIGN_PACBIO_CCS", "CONFIRMED"
    elif original_processed_fastq: classification, action, confidence = "PACBIO_PROCESSED_FASTQ", "ALIGN_PACBIO_PROCESSED", "CONFIRMED"
    elif header == "PACBIO_SUBREAD": classification, action, confidence = "PACBIO_SUBREAD_FASTQ_ONLY", "QUARANTINE", "CONFIRMED"
    elif header == "PACBIO_CCS": classification, action, confidence = "PACBIO_CCS_FASTQ", "ALIGN_PACBIO_CCS", "CONFIRMED"
    elif pacbio_declared and fastq_artifacts and fastq_subread_hint and header == "UNKNOWN":
        reasons.append("ARCHIVE_FASTQ_HEADERS_LACK_MOLECULE_GROUPING")
        classification, action, confidence = "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE", "QUARANTINE", "PROBABLE"
    elif header == "ONT": classification, action, confidence = "ONT_FASTQ", "ALIGN_ONT", "PROBABLE"
    else: classification, action, confidence = "UNKNOWN", "QUARANTINE", "INSUFFICIENT"
    if not reasons: reasons.append("NO_SUPPORTED_ARTIFACT")
    return {"submitted_representation": submitted, "header_representation": header, "classification": classification,
            "proposed_action": action, "confidence": confidence, "reason_codes": ";".join(reasons),
            "status": "QUARANTINED" if classification in ("CONFLICT", "UNKNOWN", "MIXED", "PACBIO_SUBREAD_FASTQ_ONLY") else "READY_FOR_REVIEW"}


def validate_approved(rows: list[dict]) -> None:
    allowed = {"PACBIO_CCS_FASTQ": "ALIGN_PACBIO_CCS", "PACBIO_PROCESSED_FASTQ": "ALIGN_PACBIO_PROCESSED", "PACBIO_CCS_BAM": "CANONICALISE_CCS_THEN_ALIGN", "ONT_FASTQ": "ALIGN_ONT", "PACBIO_SUBREAD_BAM": "RUN_CCS_THEN_ALIGN"}
    for n, row in enumerate(rows, 2):
        if row.get("status", "APPROVED") != "APPROVED" or row.get("review_decision") != "APPROVE": raise ValueError(f"row {n}: run is not approved")
        if row.get("classification") not in allowed or row.get("proposed_action") != allowed[row["classification"]]: raise ValueError(f"row {n}: incompatible classification/action")
        if not RUN.fullmatch(row.get("run_accession", "")) or not MD5.fullmatch(row.get("selected_artifact_md5", "")): raise ValueError(f"row {n}: unsafe accession or checksum")
        if not SAFE_URI.fullmatch(row.get("selected_artifact_uri", "")) or not SAFE_FILE.fullmatch(row.get("selected_artifact_basename", "")): raise ValueError(f"row {n}: unsafe artifact")
        if row.get("minimap2_preset") not in ("splice", "splice:hq"): raise ValueError(f"row {n}: unsupported minimap2 preset")
        if not SHA256.fullmatch(row.get("classification_report_sha256", "")): raise ValueError(f"row {n}: missing classification report hash")
        if not row.get("reviewer") or not row.get("reviewed_at"): raise ValueError(f"row {n}: missing reviewer audit fields")
        expected = row.get("expected_header_representation")
        allowed_expected = {"ONT"} if row["classification"] == "ONT_FASTQ" else ({"PACBIO_CCS", "PACBIO_CCS_ORIGINAL"} if row["classification"] == "PACBIO_CCS_FASTQ" else {"PACBIO_PROCESSED"} if row["classification"] == "PACBIO_PROCESSED_FASTQ" else {"PACBIO_CCS"})
        if expected not in allowed_expected: raise ValueError(f"row {n}: incompatible expected header representation")
        basename = row["selected_artifact_basename"].lower()
        if row["classification"] == "ONT_FASTQ" and not basename.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")): raise ValueError(f"row {n}: ONT route requires FASTQ")
        if row["classification"] == "PACBIO_CCS_FASTQ" and not basename.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")): raise ValueError(f"row {n}: CCS FASTQ route requires FASTQ")
        if row["classification"] == "PACBIO_PROCESSED_FASTQ" and not basename.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")): raise ValueError(f"row {n}: processed FASTQ route requires FASTQ")
        if row["classification"] == "PACBIO_SUBREAD_BAM" and not basename.endswith(".subreads.bam"): raise ValueError(f"row {n}: subread route requires .subreads.bam")
        if row["classification"] == "PACBIO_CCS_BAM" and not basename.endswith((".ccs.bam", ".hifi_reads.bam")): raise ValueError(f"row {n}: CCS BAM route requires a consensus BAM")
        side_uris = [value for value in row.get("required_artifact_uris", "").split(";") if value]
        side_md5s = [value for value in row.get("required_artifact_md5s", "").split(";") if value]
        if len(side_uris) != len(side_md5s): raise ValueError(f"row {n}: sidecar URI/checksum mismatch")
        for uri, checksum in zip(side_uris, side_md5s):
            if not SAFE_URI.fullmatch(uri) or not MD5.fullmatch(checksum) or not SAFE_FILE.fullmatch(Path(urllib.parse.urlparse(uri).path).name): raise ValueError(f"row {n}: unsafe sidecar")


def _cli() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("probe"); p.add_argument("fastq"); p.add_argument("output")
    f = sub.add_parser("validate-fastq"); f.add_argument("fastq"); f.add_argument("expected"); f.add_argument("output")
    v = sub.add_parser("validate-approved"); v.add_argument("manifest")
    b = sub.add_parser("build-approved"); b.add_argument("classification_tsv"); b.add_argument("report_file"); b.add_argument("output")
    i = sub.add_parser("inspect")
    i.add_argument("normalised_manifest"); i.add_argument("metadata_json"); i.add_argument("output_dir")
    args = parser.parse_args()
    if args.command == "probe":
        Path(args.output).write_text(json.dumps(probe_fastq(args.fastq), indent=2) + "\n")
    elif args.command == "validate-fastq":
        result = validate_fastq(args.fastq, args.expected)
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    elif args.command == "validate-approved":
        with open(args.manifest, newline="") as handle: rows = list(csv.DictReader(handle, delimiter="\t"))
        validate_approved(rows)
    elif args.command == "build-approved":
        report_hash = hashlib.sha256(Path(args.report_file).read_bytes()).hexdigest()
        with open(args.classification_tsv, newline="") as handle: source_rows = list(csv.DictReader(handle, delimiter="\t"))
        output_rows = []
        for row in source_rows:
            if row.get("status") != "APPROVED" or row.get("review_decision") != "APPROVE": continue
            preset = "splice" if row.get("classification") == "ONT_FASTQ" else "splice:hq"
            output_rows.append({"run_accession": row["run_accession"], "tissue": row.get("tissue", "unknown"), "description": row.get("description", "unknown"),
                                "classification": row["classification"], "proposed_action": row["proposed_action"],
                                "selected_artifact_uri": row["selected_artifact_uri"], "selected_artifact_md5": row["selected_artifact_md5"],
                                "selected_artifact_basename": row["selected_artifact_basename"], "minimap2_preset": preset,
                                "expected_header_representation": "ONT" if row["classification"] == "ONT_FASTQ" else ("PACBIO_CCS_ORIGINAL" if "NCBI_ORIGINAL_CCS_FASTQ" in row.get("reason_codes", "") else "PACBIO_CCS" if row["classification"] == "PACBIO_CCS_FASTQ" else "PACBIO_PROCESSED"),
                                "classification_report_sha256": report_hash, "reviewer": row.get("reviewer", "unknown"), "reviewed_at": row.get("reviewed_at", "unknown"),
                                "status": "APPROVED", "review_decision": "APPROVE"})
        with open(args.output, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=APPROVED_FIELDS + ["status", "review_decision"], delimiter="\t"); writer.writeheader(); writer.writerows(output_rows)
        validate_approved(output_rows)
    else:
        out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
        metadata = json.loads(Path(args.metadata_json).read_text())
        rows, artifacts, probes = [], [], []
        with open(args.normalised_manifest, newline="") as handle:
            for source_row in csv.DictReader(handle, delimiter="\t"):
                run = source_row["run_accession"]
                md = metadata.get(run, {})
                raw_artifacts = []
                if md.get("submitted_ftp"):
                    raw_artifacts.extend(expand_artifacts(run, "ENA_SUBMITTED", md.get("submitted_ftp", ""), md.get("submitted_md5", ""), md.get("submitted_format", "")))
                if md.get("fastq_ftp"):
                    raw_artifacts.extend(expand_artifacts(run, "ENA_FASTQ", md.get("fastq_ftp", ""), md.get("fastq_md5", ""), "FASTQ"))
                if not raw_artifacts and source_row.get("url"):
                    raw_artifacts = expand_artifacts(run, source_row.get("source", "ENA"), source_row["url"], source_row["md5"], "FASTQ")
                # NCBI Original Format is evidence about what the depositor
                # supplied. Its cloud URI is retained for provenance, but the
                # ENA FASTQ remains the default acquisition artifact when the
                # original is not directly downloadable here.
                for original in md.get("ncbi_original_files", []):
                    basename = original.get("filename", "")
                    if not basename:
                        continue
                    raw_artifacts.append({"run_accession": run, "source": "NCBI_ORIGINAL",
                                          "uri": original.get("uri") or f"ncbi-original://{run}/{basename}",
                                          "basename": basename, "md5": original.get("md5", ""),
                                          "format": original.get("format", artifact_role(basename)),
                                          "artifact_role": artifact_role(basename), "selected": "false"})
                probe = {"header_representation": "UNKNOWN", "records_sampled": 0, "subread_count": 0, "ccs_count": 0, "ont_uuid_count": 0, "malformed_count": 0}
                selected = next((a for a in raw_artifacts if a["artifact_role"] == "FASTQ"), None)
                if selected:
                    try:
                        probe = probe_fastq(selected["uri"])
                    except (OSError, ValueError, TimeoutError) as exc:
                        probe["probe_failure"] = str(exc)
                probes.append({"run_accession": run, **probe})
                (out / f"{run}.header_probe.json").write_text(json.dumps(probe, indent=2) + "\n")
                result = classify(declared_platform=source_row.get("platform", ""), ena_platform=md.get("instrument_platform", ""), sra_platform=md.get("sra_platform", ""),
                                  library_strategy=md.get("library_strategy", ""), library_source=md.get("library_source", ""),
                                  library_selection=md.get("library_selection", ""), artifacts=raw_artifacts, probe=probe)
                if result["classification"] in ("PACBIO_SUBREAD_BAM", "PACBIO_CCS_BAM"):
                    selected = next(a for a in raw_artifacts if a["artifact_role"] in ("SUBREADS_BAM", "CCS_BAM"))
                if selected:
                    selected["selected"] = "true"
                artifacts.extend(raw_artifacts)
                row = {field: "unknown" for field in CLASSIFICATION_FIELDS}
                row.update({"run_accession": run, "tissue": source_row.get("tissue", "unknown"), "description": source_row.get("description", "unknown"),
                            "manifest_declared_platform": source_row.get("platform", "unknown"), "ena_platform": md.get("instrument_platform", "unknown"),
                            "ena_instrument_model": md.get("instrument_model", "unknown"), "sra_platform": md.get("sra_platform", "unavailable"),
                            "sra_spot_group": md.get("sra_spot_group", "unavailable"), "selected_artifact_uri": selected["uri"] if selected else "unknown",
                            "library_strategy": md.get("library_strategy", "unknown"), "library_source": md.get("library_source", "unknown"),
                            "library_selection": md.get("library_selection", "unknown"),
                            "ncbi_original_alias": md.get("ncbi_original_alias", "unknown"),
                            "ncbi_original_representation": md.get("ncbi_original_representation", "unknown"),
                            "selected_artifact_md5": selected["md5"] if selected else "unknown", "selected_artifact_basename": selected["basename"] if selected else "unknown",
                            "header_records_sampled": probe.get("records_sampled", 0), "header_subread_count": probe.get("subread_count", 0),
                            "header_ccs_count": probe.get("ccs_count", 0), "header_ont_uuid_count": probe.get("ont_uuid_count", 0), "header_malformed_count": probe.get("malformed_count", 0)})
                row.update(result); rows.append(row)
                row["review_decision"] = "PENDING"
                (out / f"{run}.metadata.json").write_text(json.dumps(md, indent=2) + "\n")
                raw_response = md.get("ena_raw_response")
                if raw_response and Path(raw_response).is_file(): shutil.copyfile(raw_response, out / f"{run}.ena.raw.json")
        with (out / "run_classification.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CLASSIFICATION_FIELDS, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
        with (out / "run_artifacts.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=ARTIFACT_FIELDS, delimiter="\t"); writer.writeheader(); writer.writerows(artifacts)
        with (out / "header_probe.tsv").open("w", newline="") as handle:
            fields = ["run_accession", "header_representation", "records_sampled", "subread_count", "ccs_count", "ont_uuid_count", "malformed_count", "probe_failure", "probe_method", "probe_uri", "probe_bytes", "first_tokens"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
            for probe in probes:
                writer.writerow({"run_accession": probe["run_accession"], "header_representation": probe.get("header_representation", "UNKNOWN"), "records_sampled": probe.get("records_sampled", 0), "subread_count": probe.get("subread_count", 0), "ccs_count": probe.get("ccs_count", 0), "ont_uuid_count": probe.get("ont_uuid_count", 0), "malformed_count": probe.get("malformed_count", 0), "probe_failure": probe.get("probe_failure", ""), "probe_method": probe.get("probe_method", ""), "probe_uri": probe.get("probe_uri", ""), "probe_bytes": probe.get("probe_bytes", 0), "first_tokens": ";".join(probe.get("first_tokens", []))})
        with (out / "molecule_audit.tsv").open("w") as handle:
            handle.write("run_accession\trecords_sampled\tdistinct_sampled_molecules\tstatus\n")
            for probe in probes: handle.write(f"{probe['run_accession']}\t{probe.get('records_sampled', 0)}\t{probe.get('molecule_count', 0)}\tPROBE_ONLY\n")
        summary = Counter((r["classification"], r["confidence"]) for r in rows)
        (out / "classification_summary.tsv").write_text("classification\tconfidence\truns\n" + "\n".join(f"{k[0]}\t{k[1]}\t{v}" for k, v in sorted(summary.items())) + "\n")
        reasons = Counter(code for row in rows for code in row.get("reason_codes", "").split(";") if code)
        (out / "reason_code_summary.tsv").write_text("reason_code\truns\n" + "\n".join(f"{code}\t{count}" for code, count in sorted(reasons.items())) + "\n")


if __name__ == "__main__": _cli()
