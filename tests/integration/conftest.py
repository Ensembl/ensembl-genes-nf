"""
Shared fixtures for Nextflow integration tests.

Tests run each pipeline in --stub mode against minimal fixture data.
Stub mode replaces script bodies with `touch` on all output files,
verifying pipeline structure (channels, publishDir, manifest) without
executing any real compute.

Requirements:
  - nextflow >= 23.10 in PATH or NEXTFLOW_BIN env var
  - No GPU/cluster needed; uses -profile local
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
PIPELINES_DIR = Path(__file__).parent.parent.parent / "pipelines"

NEXTFLOW_BIN = os.environ.get("NEXTFLOW_BIN", str(Path.home() / ".local/bin/nextflow"))

# Java 17+ required by Nextflow; override JAVA_HOME if needed
JAVA_HOME = os.environ.get(
    "JAVA_HOME",
    "/opt/homebrew/Cellar/openjdk@21/21.0.10/libexec/openjdk.jdk/Contents/Home",
)


def run_nextflow_stub(pipeline_name: str, params: dict, tmp_path: Path) -> dict:
    """
    Run a Nextflow pipeline in --stub mode with -profile local,conda.

    Returns the parsed output_manifest.json if it exists, else {}.
    Raises subprocess.CalledProcessError on non-zero exit.
    """
    pipeline_dir = PIPELINES_DIR / pipeline_name
    outdir = tmp_path / "output"
    outdir.mkdir()

    # Build --params-file JSON to avoid shell quoting issues
    params_with_outdir = {"outdir": str(outdir), **params}
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(params_with_outdir))

    cmd = [
        NEXTFLOW_BIN,
        "run", str(pipeline_dir / "main.nf"),
        "-stub",
        "-profile", "local",
        "-params-file", str(params_file),
        "-work-dir", str(tmp_path / "work"),
        "-ansi-log", "false",
    ]

    env = os.environ.copy()
    if JAVA_HOME and Path(JAVA_HOME).exists():
        env["JAVA_HOME"] = JAVA_HOME
        env["PATH"] = str(Path(JAVA_HOME) / "bin") + ":" + env.get("PATH", "")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        env=env,
    )

    if result.returncode != 0:
        print("\n--- Nextflow stdout ---", file=sys.stderr)
        print(result.stdout[-3000:], file=sys.stderr)
        print("--- Nextflow stderr ---", file=sys.stderr)
        print(result.stderr[-3000:], file=sys.stderr)
        result.check_returncode()

    manifest_path = outdir / "output_manifest.json"
    if manifest_path.exists():
        return json.loads(manifest_path.read_text())
    return {}


@pytest.fixture(scope="session")
def nextflow_available() -> bool:
    """Skip tests if nextflow binary is not found."""
    if not Path(NEXTFLOW_BIN).exists():
        pytest.skip(f"nextflow not found at {NEXTFLOW_BIN}; set NEXTFLOW_BIN env var")
    return True


@pytest.fixture()
def minimal_gff3(tmp_path) -> Path:
    """Copy minimal.gff3 fixture into tmp_path and return its path."""
    src = FIXTURES_DIR / "minimal.gff3"
    dst = tmp_path / "minimal.gff3"
    shutil.copy(src, dst)
    return dst


@pytest.fixture()
def minimal_donor_gff3(tmp_path) -> Path:
    src = FIXTURES_DIR / "minimal_donor.gff3"
    dst = tmp_path / "minimal_donor.gff3"
    shutil.copy(src, dst)
    return dst


@pytest.fixture()
def minimal_repeats_gff3(tmp_path) -> Path:
    src = FIXTURES_DIR / "minimal_repeats.gff3"
    dst = tmp_path / "minimal_repeats.gff3"
    shutil.copy(src, dst)
    return dst


@pytest.fixture()
def minimal_fasta(tmp_path) -> Path:
    src = FIXTURES_DIR / "minimal.fa"
    dst = tmp_path / "minimal.fa"
    shutil.copy(src, dst)
    return dst
