#!/usr/bin/env python3
"""Tests for GDM inventory-driven manifest generation."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "generate_manifest.py"
SPEC = importlib.util.spec_from_file_location("generate_manifest", SCRIPT)
assert SPEC and SPEC.loader
generate_manifest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(generate_manifest)


class GenerateManifestTests(unittest.TestCase):
    def test_gdm_readme_rows_are_split_by_assembly_and_allow_scallop(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            species = root / "test_species_gca000000001v1"
            inventory = species / "readme.csv"
            bam = species / "10001" / "SRR00000001" / "alignment" / "SRR00000001_Aligned.sortedByCoord.out.bam"
            sj = bam.parent / "SRR00000001_SJ.out.tab"
            bam.parent.mkdir(parents=True)
            bam.write_text("bam\n", encoding="utf-8")
            sj.write_text("sj\n", encoding="utf-8")
            inventory.write_text(
                "taxon_id,assembly_accession,run_accession,bam_file,splice_junction_file\n"
                f"10001,GCA_000000001.1,SRR00000001,{bam.relative_to(species)},{sj.relative_to(species)}\n",
                encoding="utf-8",
            )
            output = root / "manifests"
            args = generate_manifest.parse_args([
                "--mode", "gdm", "--base-dir", str(root),
                "--gene-model-source", "scallop", "--track-types", "gene_model,splice_junction",
                "--output-dir", str(output), "--no-strict",
            ])
            rows, audit, exceptions = generate_manifest.discover(args)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["gca_accession"], "GCA_000000001.1")
            self.assertEqual(rows[0]["gtf"], "")
            self.assertEqual(rows[0]["track_types"], "gene_model,splice_junction")
            self.assertEqual(len(audit), 1)
            self.assertEqual(exceptions, [])
            generate_manifest.write_split_outputs(args, rows, audit, exceptions)
            manifest = output / "GCA_000000001.1" / "run_manifest.tsv"
            self.assertTrue(manifest.is_file())
            self.assertIn("SRR00000001", manifest.read_text(encoding="utf-8"))
            self.assertTrue((output / "gdm_manifest_index.tsv").is_file())

    def test_duplicate_inventory_rows_with_same_bam_are_collapsed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            species = root / "test_species_gca000000001v1"
            bam = species / "10001" / "SRR00000001" / "alignment" / "SRR00000001_Aligned.sortedByCoord.out.bam"
            sj = bam.parent / "SRR00000001_SJ.out.tab"
            bam.parent.mkdir(parents=True)
            bam.write_text("bam\n", encoding="utf-8")
            sj.write_text("sj\n", encoding="utf-8")
            with (species / "readme.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["taxon_id", "assembly_accession", "run_accession", "bam_file", "splice_junction_file", "tissue"])
                row = ["10001", "GCA_000000001.1", "SRR00000001", str(bam.relative_to(species)), str(sj.relative_to(species)), "brain"]
                writer.writerow(row)
                row[-1] = "liver"
                writer.writerow(row)
            args = generate_manifest.parse_args([
                "--mode", "gdm", "--base-dir", str(root), "--output-manifest", str(root / "manifest.tsv"),
            ])
            rows, audit, exceptions = generate_manifest.discover(args)
            self.assertEqual(len(rows), 1)
            self.assertIn("collapsed 2 duplicate inventory rows", audit[0]["warning"])
            self.assertEqual(exceptions, [])


if __name__ == "__main__":
    unittest.main()
