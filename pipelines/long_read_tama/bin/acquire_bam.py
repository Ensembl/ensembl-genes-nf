#!/usr/bin/env python3
import argparse
import hashlib
import os
import shutil
import time
import urllib.request
import urllib.error
from pathlib import Path


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def link_cached_file(cache, output):
    """Expose a validated cache file without copying its contents."""
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    destination.symlink_to(Path(cache).resolve())


def open_with_retries(request, attempts=6):
    """Retry transient public-object-store failures before giving up."""
    for attempt in range(attempts):
        try:
            return urllib.request.urlopen(request, timeout=120)
        except urllib.error.HTTPError as exc:
            transient = exc.code == 429 or 500 <= exc.code < 600
            if not transient or attempt == attempts - 1:
                raise
            time.sleep(min(60, 5 * (2 ** attempt)))
        except urllib.error.URLError:
            if attempt == attempts - 1:
                raise
            time.sleep(min(60, 5 * (2 ** attempt)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("expected_md5"); p.add_argument("filename"); p.add_argument("source_uri")
    p.add_argument("sidecar_uris"); p.add_argument("sidecar_md5s"); p.add_argument("cache_dir")
    p.add_argument("classification"); p.add_argument("accession"); p.add_argument("output_dir"); p.add_argument("report")
    a = p.parse_args()
    checksum_pending = a.expected_md5 == "PENDING"
    root = Path(a.cache_dir) / a.classification / ("unverified" if checksum_pending else a.expected_md5)
    root.mkdir(parents=True, exist_ok=True); cache = root / a.filename; lock = Path(str(cache) + ".lock")
    lock_acquired = False
    if not cache.exists():
        while True:
            try:
                lock.mkdir()
                lock_acquired = True
                break
            except FileExistsError:
                time.sleep(1)
    try:
        cache_hit = cache.exists()
        if not cache_hit:
            part = Path(str(cache) + ".part")
            # Some public object stores return 403 to Python's default
            # User-Agent from HPC task containers.  Keep the acquisition
            # self-contained: an explicit User-Agent is enough for PacBio,
            # and the Range request resumes a partial task-local transfer.
            offset = part.stat().st_size if part.exists() else 0
            headers = {"User-Agent": "Mozilla/5.0 (Ensembl long-read pipeline)"}
            if offset:
                headers["Range"] = f"bytes={offset}-"
            request = urllib.request.Request(a.source_uri, headers=headers)
            with open_with_retries(request) as src:
                resumed = offset > 0 and getattr(src, "status", None) == 206
                if offset and not resumed:
                    part.unlink()
                with open(part, "ab" if resumed else "wb") as dst:
                    shutil.copyfileobj(src, dst, length=1024 * 1024)
            actual = md5(part)
            if not checksum_pending and actual != a.expected_md5: raise SystemExit(f"BAM checksum mismatch for {a.accession}")
            os.replace(part, cache)
        # Known-checksum cache paths are content-addressed.  Do not re-read
        # multi-GB cache hits before linking them into the task work directory.
        # For PENDING inputs no checksum is available by definition.
        actual = "PENDING" if checksum_pending else a.expected_md5
        out = Path(a.output_dir); out.mkdir(exist_ok=True); link_cached_file(cache, out / a.filename)
        uris = [x for x in a.sidecar_uris.split(";") if x]
        checks = [x for x in a.sidecar_md5s.split(";") if x]
        if len(uris) != len(checks): raise SystemExit("sidecar URI/checksum mismatch")
        for uri, expected in zip(uris, checks):
            name = Path(uri.split("?", 1)[0]).name; side = Path(a.cache_dir) / a.classification / expected / name
            side.parent.mkdir(parents=True, exist_ok=True)
            side_hit = side.exists()
            if not side_hit:
                side_request = urllib.request.Request(uri, headers={"User-Agent": "Mozilla/5.0 (Ensembl long-read pipeline)"})
                with open_with_retries(side_request) as src, open(str(side) + ".part", "wb") as dst: shutil.copyfileobj(src, dst)
                if md5(str(side) + ".part") != expected: raise SystemExit(f"Sidecar checksum mismatch for {a.accession}: {name}")
                os.replace(str(side) + ".part", side)
            link_cached_file(side, out / name)
        Path(a.report).write_text("run_accession\texpected_md5\tactual_md5\tfilename\n" + f"{a.accession}\t{a.expected_md5}\t{actual}\t{a.filename}\n")
    finally:
        if lock_acquired and lock.exists():
            lock.rmdir()


if __name__ == "__main__": main()
