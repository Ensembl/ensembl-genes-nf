#!/usr/bin/env python3
import argparse
import hashlib
import os
import shutil
import subprocess
import time
from pathlib import Path


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("expected_md5"); p.add_argument("filename"); p.add_argument("source_uri")
    p.add_argument("sidecar_uris"); p.add_argument("sidecar_md5s"); p.add_argument("cache_dir")
    p.add_argument("classification"); p.add_argument("accession"); p.add_argument("output_dir"); p.add_argument("report")
    a = p.parse_args()
    checksum_pending = a.expected_md5 == "PENDING"
    root = Path(a.cache_dir) / a.classification / ("unverified" if checksum_pending else a.expected_md5)
    root.mkdir(parents=True, exist_ok=True); cache = root / a.filename; lock = Path(str(cache) + ".lock")
    while True:
        try: lock.mkdir(); break
        except FileExistsError: time.sleep(1)
    try:
        if not cache.exists():
            part = Path(str(cache) + ".part")
            # Some public object stores return 403 to Python urllib from HPC
            # task containers while accepting curl.  Use curl for resumable,
            # retried transfer and keep the checksum gate below authoritative.
            subprocess.run([
                "curl", "--fail", "--location", "--retry", "5", "--retry-all-errors",
                "--retry-delay", "5", "--connect-timeout", "60", "--continue-at", "-",
                "--output", str(part), a.source_uri,
            ], check=True)
            actual = md5(part)
            if not checksum_pending and actual != a.expected_md5: raise SystemExit(f"BAM checksum mismatch for {a.accession}")
            os.replace(part, cache)
        actual = md5(cache)
        if not checksum_pending and actual != a.expected_md5: raise SystemExit(f"BAM checksum mismatch for {a.accession}")
        out = Path(a.output_dir); out.mkdir(exist_ok=True); shutil.copy2(cache, out / a.filename)
        uris = [x for x in a.sidecar_uris.split(";") if x]
        checks = [x for x in a.sidecar_md5s.split(";") if x]
        if len(uris) != len(checks): raise SystemExit("sidecar URI/checksum mismatch")
        for uri, expected in zip(uris, checks):
            name = Path(uri.split("?", 1)[0]).name; side = Path(a.cache_dir) / a.classification / expected / name
            side.parent.mkdir(parents=True, exist_ok=True)
            if not side.exists():
                with urlopen(uri) as src, open(str(side) + ".part", "wb") as dst: shutil.copyfileobj(src, dst)
                if md5(str(side) + ".part") != expected: raise SystemExit(f"Sidecar checksum mismatch for {a.accession}: {name}")
                os.replace(str(side) + ".part", side)
            if md5(side) != expected: raise SystemExit(f"Cached sidecar checksum mismatch for {a.accession}: {name}")
            shutil.copy2(side, out / name)
        Path(a.report).write_text("run_accession\texpected_md5\tactual_md5\tfilename\n" + f"{a.accession}\t{a.expected_md5}\t{actual}\t{a.filename}\n")
    finally: lock.rmdir()


if __name__ == "__main__": main()
