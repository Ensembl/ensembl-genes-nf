"""
Stub-mode integration tests for Nextflow pipelines.

Each test runs a pipeline with --stub, which executes the `stub:` block
instead of the `script:` block. This validates:

  1. The pipeline compiles and its process DAG is wired correctly
  2. All expected output files land in their publishDir
  3. output_manifest.json is valid JSON with the right pipeline name
  4. Channel connections between processes are correct (stub GFF3 → manifest input)

These tests do NOT validate biological correctness — they verify pipeline
structure and plumbing only.

Run with:
  cd ensembl-genes-nf
  pytest tests/integration/ -v --tb=short
  NEXTFLOW_BIN=/path/to/nextflow pytest tests/integration/ -v
"""

import json
import shutil
from pathlib import Path

import pytest

from conftest import FIXTURES_DIR, PIPELINES_DIR, run_nextflow_stub


# ---------------------------------------------------------------------------
# utr_addition
# ---------------------------------------------------------------------------

class TestUtrAdditionStub:
    """Stub-mode tests for the utr_addition pipeline."""

    def _run(self, minimal_gff3, minimal_donor_gff3, tmp_path):
        return run_nextflow_stub(
            "utr_addition",
            params={
                "consolidated_gff3": str(minimal_gff3),
                "donor_gff3_files":  str(minimal_donor_gff3),
            },
            tmp_path=tmp_path,
        )

    def test_pipeline_completes_in_stub_mode(self, nextflow_available, minimal_gff3, minimal_donor_gff3, tmp_path):
        """Pipeline exits 0 in stub mode."""
        manifest = self._run(minimal_gff3, minimal_donor_gff3, tmp_path)
        assert isinstance(manifest, dict)

    def test_manifest_is_valid_json(self, nextflow_available, minimal_gff3, minimal_donor_gff3, tmp_path):
        """output_manifest.json exists and is parseable JSON."""
        manifest_path = tmp_path / "output" / "output_manifest.json"
        self._run(minimal_gff3, minimal_donor_gff3, tmp_path)
        assert manifest_path.exists(), "output_manifest.json not found"
        data = json.loads(manifest_path.read_text())
        assert data.get("pipeline") == "utr_addition"

    def test_gff3_published_to_utr_addition_subdir(self, nextflow_available, minimal_gff3, minimal_donor_gff3, tmp_path):
        """ADD_UTRS publishes a GFF3 stub to outdir/utr_addition/."""
        self._run(minimal_gff3, minimal_donor_gff3, tmp_path)
        utr_dir = tmp_path / "output" / "utr_addition"
        assert utr_dir.exists(), f"publishDir utr_addition/ not found under {tmp_path / 'output'}"
        gff3_files = list(utr_dir.glob("*.gff3"))
        assert len(gff3_files) >= 1, f"No GFF3 files in utr_addition/: {list(utr_dir.iterdir())}"

    def test_stub_gff3_contains_gff3_version_header(self, nextflow_available, minimal_gff3, minimal_donor_gff3, tmp_path):
        """The stub GFF3 starts with ##gff-version 3 (written by stub block)."""
        self._run(minimal_gff3, minimal_donor_gff3, tmp_path)
        utr_dir = tmp_path / "output" / "utr_addition"
        gff3 = next(utr_dir.glob("*.gff3"))
        content = gff3.read_text()
        assert content.startswith("##gff-version 3"), f"GFF3 header missing, got: {content[:50]!r}"

    def test_stub_gff3_has_utr_addition_source(self, nextflow_available, minimal_gff3, minimal_donor_gff3, tmp_path):
        """The stub GFF3 contains a gene feature (written by ADD_UTRS stub)."""
        self._run(minimal_gff3, minimal_donor_gff3, tmp_path)
        utr_dir = tmp_path / "output" / "utr_addition"
        gff3 = next(utr_dir.glob("*.gff3"))
        content = gff3.read_text()
        assert "gene" in content, f"No gene feature in stub GFF3: {content!r}"

    def test_output_filename_derived_from_input(self, nextflow_available, minimal_gff3, minimal_donor_gff3, tmp_path):
        """Output filename is <input_stem>.with_utrs.gff3."""
        self._run(minimal_gff3, minimal_donor_gff3, tmp_path)
        utr_dir = tmp_path / "output" / "utr_addition"
        gff3 = next(utr_dir.glob("*.gff3"))
        assert gff3.name.endswith(".with_utrs.gff3"), f"Unexpected output name: {gff3.name}"


# ---------------------------------------------------------------------------
# finalise_geneset
# ---------------------------------------------------------------------------

class TestFinaliseGenesetStub:
    """Stub-mode tests for the finalise_geneset pipeline."""

    def _run(self, minimal_gff3, minimal_repeats_gff3, tmp_path, extra=None):
        params = {
            "input_gff3":  str(minimal_gff3),
            "repeat_gff3": str(minimal_repeats_gff3),
        }
        if extra:
            params.update(extra)
        return run_nextflow_stub("finalise_geneset", params=params, tmp_path=tmp_path)

    def test_pipeline_completes_in_stub_mode(self, nextflow_available, minimal_gff3, minimal_repeats_gff3, tmp_path):
        self._run(minimal_gff3, minimal_repeats_gff3, tmp_path)

    def test_manifest_is_valid_json(self, nextflow_available, minimal_gff3, minimal_repeats_gff3, tmp_path):
        manifest_path = tmp_path / "output" / "output_manifest.json"
        self._run(minimal_gff3, minimal_repeats_gff3, tmp_path)
        assert manifest_path.exists()
        data = json.loads(manifest_path.read_text())
        assert data.get("pipeline") == "finalise_geneset"

    def test_gff3_published_to_finalise_geneset_subdir(self, nextflow_available, minimal_gff3, minimal_repeats_gff3, tmp_path):
        """SELECT_CANONICAL publishes GFF3 to outdir/finalise_geneset/."""
        self._run(minimal_gff3, minimal_repeats_gff3, tmp_path)
        fg_dir = tmp_path / "output" / "finalise_geneset"
        assert fg_dir.exists(), "finalise_geneset publishDir not found"
        gff3_files = list(fg_dir.glob("*.gff3"))
        assert len(gff3_files) >= 1, f"No GFF3 in finalise_geneset/: {list(fg_dir.iterdir())}"

    def test_pipeline_handles_absent_selenoprotein_fasta(self, nextflow_available, minimal_gff3, minimal_repeats_gff3, tmp_path):
        """Pipeline runs without selenoprotein_fasta (uses NO_FILE default)."""
        self._run(minimal_gff3, minimal_repeats_gff3, tmp_path)  # no selenoprotein_fasta

    def test_all_five_processes_produce_outputs(self, nextflow_available, minimal_gff3, minimal_repeats_gff3, tmp_path):
        """
        All 5 processes in the pipeline chain (filter → pseudogenes → readthrough
        → selenoproteins → canonical) complete, producing a single final GFF3.
        """
        self._run(minimal_gff3, minimal_repeats_gff3, tmp_path)
        fg_dir = tmp_path / "output" / "finalise_geneset"
        # The canonical selector is the last process; only its output survives
        gff3_files = list(fg_dir.glob("*.gff3"))
        assert len(gff3_files) >= 1


# ---------------------------------------------------------------------------
# consolidate
# ---------------------------------------------------------------------------

class TestConsolidateStub:
    """Stub-mode tests for the consolidate pipeline."""

    def _gff3_file(self, tmp_path) -> Path:
        """Copy minimal.gff3 fixture into tmp_path and return its path."""
        dst = tmp_path / "rnaseq.merged.gff3"
        shutil.copy(FIXTURES_DIR / "minimal.gff3", dst)
        return dst

    def _params(self, gff3_file):
        # Use gff3_files (explicit path list) rather than gff3_dir to avoid
        # the **/*.gff3 recursive glob which requires at least one subdir level.
        return {"gff3_files": str(gff3_file)}

    def test_pipeline_completes_in_stub_mode(self, nextflow_available, tmp_path):
        gff3_file = self._gff3_file(tmp_path)
        run_nextflow_stub("consolidate", params=self._params(gff3_file), tmp_path=tmp_path)

    def test_manifest_is_valid_json(self, nextflow_available, tmp_path):
        gff3_file = self._gff3_file(tmp_path)
        manifest_path = tmp_path / "output" / "output_manifest.json"
        run_nextflow_stub("consolidate", params=self._params(gff3_file), tmp_path=tmp_path)
        assert manifest_path.exists()
        data = json.loads(manifest_path.read_text())
        assert data.get("pipeline") == "consolidate"

    def test_gff3_published_to_consolidate_subdir(self, nextflow_available, tmp_path):
        """Consolidated GFF3 lands in outdir/consolidate/."""
        gff3_file = self._gff3_file(tmp_path)
        run_nextflow_stub("consolidate", params=self._params(gff3_file), tmp_path=tmp_path)
        con_dir = tmp_path / "output" / "consolidate"
        assert con_dir.exists(), "consolidate publishDir not found"
        gff3_files = list(con_dir.glob("*.gff3"))
        assert len(gff3_files) >= 1


# ---------------------------------------------------------------------------
# Pipeline chain integrity: utr_addition → finalise_geneset
# ---------------------------------------------------------------------------

class TestPipelineChainIntegrity:
    """
    Verify that the two post-consolidation pipelines can be chained:
    the output filename convention from utr_addition must be consumable
    as input to finalise_geneset.
    """

    def test_utr_addition_output_can_feed_finalise_geneset(
        self, nextflow_available, minimal_gff3, minimal_donor_gff3, minimal_repeats_gff3, tmp_path
    ):
        """
        Run utr_addition in stub mode, take its output GFF3, and feed it
        directly into finalise_geneset stub mode. Both should succeed.
        """
        # Stage 1: UTR addition
        utr_tmp = tmp_path / "utr"
        utr_tmp.mkdir()
        run_nextflow_stub(
            "utr_addition",
            params={
                "consolidated_gff3": str(minimal_gff3),
                "donor_gff3_files":  str(minimal_donor_gff3),
            },
            tmp_path=utr_tmp,
        )
        utr_dir = utr_tmp / "output" / "utr_addition"
        utr_gff3 = next(utr_dir.glob("*.gff3"))

        # Stage 2: Finalise geneset using the UTR-added GFF3
        fg_tmp = tmp_path / "fg"
        fg_tmp.mkdir()
        run_nextflow_stub(
            "finalise_geneset",
            params={
                "input_gff3":  str(utr_gff3),
                "repeat_gff3": str(minimal_repeats_gff3),
            },
            tmp_path=fg_tmp,
        )

        fg_dir = fg_tmp / "output" / "finalise_geneset"
        assert fg_dir.exists()
        final_gff3_files = list(fg_dir.glob("*.gff3"))
        assert len(final_gff3_files) >= 1, "finalise_geneset produced no output"
