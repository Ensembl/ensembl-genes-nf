#!/usr/bin/env python3
"""Run executable smoke checks for every pinned long-read container."""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shlex
import signal
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


COMMANDS = {
    "minimap2": "minimap2 --version",
    "samtools": "samtools --version",
    "tama": "tama_collapse.py -v 1",
    "stringtie2": "stringtie --version",
    "stringtie3": "stringtie --version",
    "tmerge": "perl -e 'print \"$^V\\n\"'",
    "pbccs": "ccs --version",
    "isoquant": "isoquant.py --version",
    "flair": "flair --version",
    "bambu": "Rscript -e 'cat(as.character(packageVersion(\"bambu\")), \"\\n\")'",
    "transgram": "/TransGram/TransGram --version 2>&1 || test $? -eq 1",
    "diamond": "diamond version",
    "python_helper": "python3 --version",
}
REF_RE = re.compile(r"^([A-Za-z0-9_.-]+):\s+([^\s#]+)")


def text_value(value: str | bytes | None) -> str:
    """Normalize subprocess output across Python versions and timeout paths."""
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value


def run_command(argv: list[str], timeout: int) -> tuple[int, str, str]:
    """Run one smoke command and terminate descendants when it times out."""
    completed = subprocess.Popen(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=(os.name == "posix"),
    )
    try:
        stdout, stderr = completed.communicate(timeout=timeout)
        return completed.returncode, text_value(stdout), text_value(stderr)
    except subprocess.TimeoutExpired as error:
        try:
            if os.name == "posix":
                os.killpg(completed.pid, signal.SIGKILL)
            else:
                completed.kill()
        except ProcessLookupError:
            pass
        try:
            completed.wait(timeout=5)
        except subprocess.TimeoutExpired:
            completed.kill()
            completed.wait()
        if completed.stdout is not None:
            completed.stdout.close()
        if completed.stderr is not None:
            completed.stderr.close()
        cleanup_timed_out_container(argv)
        return (
            124,
            text_value(error.stdout),
            text_value(error.stderr) + f"\nTimed out after {timeout}s\n",
        )


def cleanup_timed_out_container(argv: list[str]) -> None:
    """Best-effort cleanup for Apptainer commands launched through Lima.

    ``limactl shell`` can leave the guest-side Apptainer process alive after
    the host SSH wrapper is killed.  Only target the exact immutable image
    reference from this smoke command, and only when the configured engine is
    a wrapper that explicitly names a Lima instance.
    """
    if not argv:
        return
    engine = Path(shutil.which(argv[0]) or argv[0])
    try:
        wrapper = engine.read_text()
    except (OSError, UnicodeError):
        return
    instance_match = re.search(r"limactl\s+shell\s+([A-Za-z0-9_.-]+)", wrapper)
    reference = next((item for item in argv if item.startswith("docker://")), None)
    if not instance_match or not reference:
        return
    cleanup = [
        "limactl",
        "shell",
        instance_match.group(1),
        "--",
        "sh",
        "-lc",
        f"pkill -KILL -f {shlex.quote(reference)} 2>/dev/null || true",
    ]
    try:
        subprocess.run(cleanup, check=False, timeout=20, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.TimeoutExpired):
        pass


def read_registry(path: Path) -> dict[str, str]:
    refs = {}
    for line in path.read_text().splitlines():
        match = REF_RE.match(line.strip())
        if match:
            refs[match.group(1)] = match.group(2)
    missing = sorted(set(COMMANDS) - set(refs))
    if missing:
        raise ValueError(f"container registry is missing smoke targets: {', '.join(missing)}")
    return refs


def build_argv(
    engine: str, reference: str, command: str, platform_name: str, disable_cache: bool
) -> list[str]:
    """Build an engine-specific executable/version command."""
    engine_name = Path(engine).name
    if engine_name in {"docker", "podman", "nerdctl"}:
        argv = [engine, "run", "--rm", "--platform", platform_name]
        if disable_cache:
            argv.extend(["--pull", "always"])
        argv.extend([reference, "sh", "-lc", command])
        return argv

    argv = [engine, "exec", "--cleanenv"]
    if disable_cache:
        argv.append("--disable-cache")
    image = reference if reference.startswith(("http://", "https://")) else f"docker://{reference}"
    argv.extend([image, "sh", "-lc", command])
    return argv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("registry", type=Path)
    parser.add_argument("--engine", default="apptainer")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    parser.add_argument("--platform", default="linux/amd64")
    parser.add_argument(
        "--disable-cache",
        action="store_true",
        help="force a fresh Apptainer conversion or Docker pull instead of using cache",
    )
    args = parser.parse_args()

    refs = read_registry(args.registry)
    results = []
    for name in COMMANDS:
        command = COMMANDS[name]
        image = (
            refs[name]
            if Path(args.engine).name in {"docker", "podman", "nerdctl"}
            else f"docker://{refs[name]}"
        )
        argv = build_argv(args.engine, refs[name], command, args.platform, args.disable_cache)
        started = datetime.now(timezone.utc).isoformat()
        try:
            returncode, stdout, stderr = run_command(argv, args.timeout_seconds)
        except OSError as error:
            returncode = 127
            stdout = ""
            stderr = f"Unable to execute {args.engine}: {error}\n"
        results.append(
            {
                "name": name,
                "reference": refs[name],
                "image": image,
                "platform": args.platform,
                "command": command,
                "argv": shlex.join(argv),
                "started_at": started,
                "returncode": returncode,
                "stdout": stdout,
                "stderr": stderr,
            }
        )

    payload = {
        "engine": args.engine,
        "platform": args.platform,
        "host": {"system": platform.system(), "machine": platform.machine()},
        "registry": str(args.registry),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0 if all(item["returncode"] == 0 for item in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
