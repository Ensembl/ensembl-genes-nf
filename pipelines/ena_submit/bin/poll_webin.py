#!/usr/bin/env python3
import argparse, base64, json, os, sys, time, urllib.request, urllib.error, re
from pathlib import Path


def basic_auth_header(username: str, password: str) -> str:
    token = base64.b64encode(f"{username}:{password}".encode()).decode()
    return f"Basic {token}"


def http_get(url: str, auth_header: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"Authorization": auth_header})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def parse_queue_json(path: Path):
    if not path.exists():
        raise RuntimeError(f"Queue response does not exist: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(
            f"Queue response is empty: {path}; the submit step likely failed or returned no body"
        )
    try:
        with path.open() as fh:
            data = json.load(fh)
    except Exception as e:
        raise RuntimeError(f"Failed to parse JSON from {path}: {e}")
    poll_url = (
        data.get("_links", {}).get("poll", {}).get("href")
        or data.get("href")
        or ""
    )
    sub_id = data.get("submissionId") or ""
    if not poll_url:
        raise RuntimeError(f"No poll href in {path}")
    return sub_id, poll_url


def extract_accession(receipt_xml: str) -> str:
    m = re.search(r'accession="([^"]+)"', receipt_xml)
    return m.group(1) if m else ""


def looks_idempotent_exists(text: str) -> bool:
    t = text.lower()
    return (
        'already exists' in t or
        'object being added already exists' in t or
        'alias not unique' in t or
        'the object you are trying to submit already exists' in t
    )


def main():
    ap = argparse.ArgumentParser(description="Poll ENA Webin queue and produce accessions.tsv")
    ap.add_argument("queues", nargs="+", help="List of *.queue.json files")
    ap.add_argument("--webin-user", required=True)
    password_group = ap.add_mutually_exclusive_group()
    password_group.add_argument("--webin-password", help="Password (discouraged: visible in task command files)")
    password_group.add_argument("--webin-password-env", default="ENA_WEBIN_PASSWORD")
    ap.add_argument("--interval", type=int, default=20)
    ap.add_argument("--max-attempts", type=int, default=30)
    ap.add_argument("--out", default="accessions.tsv")
    args = ap.parse_args()

    password = args.webin_password or os.environ.get(args.webin_password_env)
    if not password:
        ap.error(f"No Webin password supplied; set ${args.webin_password_env} or use --webin-password")
    auth = basic_auth_header(args.webin_user, password)

    failures = 0
    with open(args.out, "w") as out:
        out.write("analysis_id\tsubmission_id\taccession\tstatus\n")
        for qpath in args.queues:
            qfile = Path(qpath)
            analysis_id = qfile.stem.replace(".queue", "") if qfile.name.endswith(".queue.json") else qfile.stem
            try:
                submission_id, poll_url = parse_queue_json(qfile)
            except Exception as e:
                print(f"ERROR: {e}", file=sys.stderr)
                failures = 1
                out.write(f"{analysis_id}\t\t\tFAILED\n")
                continue

            print(f"Polling {submission_id or '?'} ({analysis_id}) ...", file=sys.stderr)
            receipt = None
            for i in range(1, args.max_attempts + 1):
                try:
                    body = http_get(poll_url, auth)
                except urllib.error.HTTPError as e:
                    # Treat 404/5xx as pending within window
                    body = e.read() or b""
                except Exception as e:
                    body = b""
                text = body.decode(errors="replace")
                if "<RECEIPT" in text:
                    receipt = text
                    break
                print(f"  attempt {i}: pending", file=sys.stderr)
                time.sleep(args.interval)

            if not receipt:
                print(f"ERROR: timed out polling {submission_id or '?'}", file=sys.stderr)
                failures = 1
                out.write(f"{analysis_id}\t{submission_id}\t\tTIMEOUT\n")
                continue

            if 'success="false"' in receipt:
                if looks_idempotent_exists(receipt):
                    acc = extract_accession(receipt)
                    print(f"EXISTS: {submission_id or '?'} (alias already registered){' -> ' + acc if acc else ''}", file=sys.stderr)
                    out.write(f"{analysis_id}\t{submission_id}\t{acc}\tEXISTS\n")
                    continue
                else:
                    print(f"ERROR: {submission_id or '?'} rejected by ENA:", file=sys.stderr)
                    print(receipt, file=sys.stderr)
                    failures = 1
                    out.write(f"{analysis_id}\t{submission_id}\t\tFAILED\n")
                    continue

            acc = extract_accession(receipt)
            print(f"OK: {submission_id or '?'} -> {acc or '?'}", file=sys.stderr)
            out.write(f"{analysis_id}\t{submission_id}\t{acc}\tSUCCESS\n")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
