#!/usr/bin/env python3
"""Fetch a taxon-scoped UniProtKB FASTA with pagination and provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


NEXT_LINK = re.compile(r"<([^>]+)>;\s*rel=\"next\"")
RETRYABLE = {429, 500, 502, 503, 504}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--api", required=True)
    parser.add_argument("--build-id", required=True)
    parser.add_argument("--source", choices=("uniprot", "swissprot"), required=True)
    parser.add_argument("--scope", choices=("species", "clade"), required=True)
    parser.add_argument("--taxon-id", required=True)
    parser.add_argument("--taxon-name", required=True)
    parser.add_argument("--query-extra", default="")
    parser.add_argument("--include-isoforms", choices=("true", "false"), default="false")
    parser.add_argument("--page-size", type=int, choices=range(1, 501), default=500)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    return parser.parse_args()


def query_for(args: argparse.Namespace) -> str:
    taxon_field = "organism_id" if args.scope == "species" else "lineage"
    query = f"({taxon_field}:{args.taxon_id})"
    if args.source == "swissprot":
        query += " AND (reviewed:true)"
    if args.query_extra:
        query += f" AND ({args.query_extra})"
    return query


def request_page(url: str):
    request = Request(
        url,
        headers={
            "Accept": "text/plain; format=fasta",
            "User-Agent": "ensembl-genes-nf/protein-db",
        },
    )
    for attempt in range(6):
        try:
            return urlopen(request, timeout=180)
        except HTTPError as exc:
            if exc.code not in RETRYABLE or attempt == 5:
                raise
            time.sleep(2**attempt)
        except URLError:
            if attempt == 5:
                raise
            time.sleep(2**attempt)
    raise RuntimeError(f"Unable to retrieve {url}")


def fetch_fasta(args: argparse.Namespace, query: str) -> tuple[int, str, str, str]:
    params = {
        "query": query,
        "format": "fasta",
        "compressed": "false",
        "includeIsoform": args.include_isoforms,
        "size": str(args.page_size),
    }
    next_url = f"{args.api}?{urlencode(params)}"
    total = "unknown"
    release = "unknown"
    sequence_count = 0
    digest = hashlib.sha256()

    with args.fasta.open("wb") as output:
        while next_url:
            with request_page(next_url) as response:
                if total == "unknown":
                    total = response.headers.get("x-total-results", "unknown")
                    release = response.headers.get("x-uniprot-release", "unknown")
                body = response.read()
                if body and not body.endswith(b"\n"):
                    body += b"\n"
                output.write(body)
                digest.update(body)
                sequence_count += body.count(b"\n>")
                if body.startswith(b">"):
                    sequence_count += 1
                match = NEXT_LINK.search(response.headers.get("Link", ""))
                next_url = match.group(1) if match else ""

    return sequence_count, total, release, digest.hexdigest()


def write_metadata(
    args: argparse.Namespace,
    query: str,
    sequence_count: int,
    total: str,
    release: str,
    digest: str,
) -> None:
    rows = [
        ("build_id", args.build_id),
        ("source", args.source),
        ("scope", args.scope),
        ("taxon_id", args.taxon_id),
        ("taxon_name", args.taxon_name),
        ("query", query),
        ("query_extra", args.query_extra),
        ("api_endpoint", args.api),
        ("uniprot_release", release),
        ("retrieved_at_utc", datetime.now(timezone.utc).isoformat()),
        ("include_isoforms", args.include_isoforms),
        ("sequence_count", str(sequence_count)),
        ("api_total_results", total),
        ("fasta_sha256", digest),
    ]
    with args.metadata.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    query = query_for(args)
    sequence_count, total, release, digest = fetch_fasta(args, query)
    if sequence_count == 0:
        raise RuntimeError(f"UniProt returned no FASTA records for query: {query}")
    write_metadata(args, query, sequence_count, total, release, digest)
    print(
        f"Fetched {sequence_count} FASTA records for {args.build_id} "
        f"(API total: {total})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
