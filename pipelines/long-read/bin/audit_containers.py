#!/usr/bin/env python3
"""Audit the long-read container registry references, including OCI digests."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ACCEPT = ", ".join(
    (
        "application/vnd.oci.image.index.v1+json",
        "application/vnd.docker.distribution.manifest.list.v2+json",
        "application/vnd.oci.image.manifest.v1+json",
        "application/vnd.docker.distribution.manifest.v2+json",
    )
)
REF_RE = re.compile(r"(?m)(?:^|['\"\s])((?:quay\.io|docker\.io)/[^'\"\s]+@sha256:[0-9a-f]{64})(?:['\"\s]|$)")


def get_json(url: str, headers: dict[str, str]) -> tuple[dict, dict[str, str]]:
    try:
        with urlopen(Request(url, headers=headers), timeout=45) as response:
            return json.load(response), dict(response.headers.items())
    except HTTPError as error:
        if error.code != 401:
            raise
        challenge = error.headers.get("WWW-Authenticate", "")
        fields = dict(re.findall(r'(\w+)="([^"]+)"', challenge))
        query = urlencode({"service": fields.get("service", ""), "scope": fields.get("scope", "")})
        with urlopen(f"{fields['realm']}?{query}", timeout=45) as response:
            token = json.load(response)["token"]
        authorized = dict(headers)
        authorized["Authorization"] = f"Bearer {token}"
        with urlopen(Request(url, headers=authorized), timeout=45) as response:
            return json.load(response), dict(response.headers.items())


def audit(ref: str) -> str:
    registry, remainder = ref.split("/", 1)
    endpoint_registry = "registry-1.docker.io" if registry == "docker.io" else registry
    repository, digest = remainder.rsplit("@", 1)
    if ":" in repository.rsplit("/", 1)[-1]:
        repository = repository.rsplit(":", 1)[0]
    endpoint = f"https://{endpoint_registry}/v2/{repository}/manifests/{digest}"
    manifest, headers = get_json(endpoint, {"Accept": ACCEPT, "User-Agent": "ensembl-long-read-container-audit"})
    if manifest.get("schemaVersion") != 2:
        raise RuntimeError("unexpected OCI/Docker manifest schema")
    if "manifests" in manifest:
        platforms = {
            ((item.get("platform") or {}).get("os"), (item.get("platform") or {}).get("architecture"))
            for item in manifest["manifests"]
        }
        if ("linux", "amd64") not in platforms:
            raise RuntimeError(f"no linux/amd64 descriptor; available={sorted(platforms)}")
        return "linux/amd64 descriptor"
    config_digest = (manifest.get("config") or {}).get("digest")
    if not config_digest:
        raise RuntimeError("manifest has no config digest")
    blob_url = f"https://{endpoint_registry}/v2/{repository}/blobs/{config_digest}"
    config, _ = get_json(blob_url, {"Accept": "application/json", **({"Authorization": headers["Authorization"]} if "Authorization" in headers else {})})
    if (config.get("os"), config.get("architecture")) != ("linux", "amd64"):
        raise RuntimeError(f"single manifest is {config.get('os')}/{config.get('architecture')}")
    return "linux/amd64 single manifest"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pipeline", type=Path)
    parser.add_argument("--platform", default="linux/amd64")
    args = parser.parse_args()
    text = "\n".join(path.read_text(errors="replace") for path in args.pipeline.rglob("*") if path.is_file() and path.suffix in {".nf", ".config", ".yml"})
    refs = sorted(set(REF_RE.findall(text)))
    if not refs:
        print("no immutable OCI references found", file=sys.stderr)
        return 2
    failed = False
    for ref in refs:
        try:
            print(f"OK\t{ref}\t{audit(ref)}")
        except Exception as error:  # noqa: BLE001 - report every image in one audit
            failed = True
            print(f"FAIL\t{ref}\t{error}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
