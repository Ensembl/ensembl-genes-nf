#!/usr/bin/env python3
import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import time
from pathlib import Path
from urllib.request import urlopen


def checksum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_gzip(path):
    with gzip.open(path, "rb") as handle:
        while handle.read(1024 * 1024):
            pass


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run"); p.add_argument("expected_md5"); p.add_argument("filename")
    p.add_argument("source_uri"); p.add_argument("cache_dir"); p.add_argument("classification"); p.add_argument("output_dir")
    p.add_argument("cpus", type=int); p.add_argument("output_checksum")
    a = p.parse_args()
    cache = Path(a.cache_dir) / a.classification / a.expected_md5 / a.filename
    cache.parent.mkdir(parents=True, exist_ok=True)
    lock = Path(str(cache) + ".lock")
    while True:
        try:
            lock.mkdir(); break
        except FileExistsError:
            time.sleep(1)
    try:
        if cache.exists():
            validate_gzip(cache)
            actual = checksum(cache)
            if actual != a.expected_md5: raise SystemExit(f"Cached FASTQ checksum mismatch for {a.run}: expected {a.expected_md5}, observed {actual}")
        else:
            acquisition = Path("acquisition"); acquisition.mkdir(exist_ok=True)
            if a.source_uri:
                fetched = acquisition / a.filename
                with urlopen(a.source_uri) as source, open(fetched, "wb") as target:
                    shutil.copyfileobj(source, target)
            else:
                subprocess.run(["fastq-dl", "-a", a.run, "--cpus", str(a.cpus)], cwd=acquisition, check=True)
                files = list(acquisition.glob("*.fastq.gz"))
                if len(files) != 1: raise SystemExit(f"Expected one long-read FASTQ for {a.run}, found {len(files)}")
                fetched = files[0]
            validate_gzip(fetched)
            if fetched.name != a.filename: raise SystemExit(f"Downloaded FASTQ filename mismatch for {a.run}")
            actual = checksum(fetched)
            if actual != a.expected_md5: raise SystemExit(f"Downloaded FASTQ checksum mismatch for {a.run}: expected {a.expected_md5}, observed {actual}")
            os.replace(fetched, cache)
        out = Path(a.output_dir); out.mkdir(exist_ok=True)
        shutil.copy2(cache, out / a.filename)
        Path(a.output_checksum).write_text("run_accession\texpected_md5\tactual_md5\tfilename\n" + f"{a.run}\t{a.expected_md5}\t{checksum(cache)}\t{a.filename}\n")
    finally:
        lock.rmdir()


if __name__ == "__main__":
    main()
