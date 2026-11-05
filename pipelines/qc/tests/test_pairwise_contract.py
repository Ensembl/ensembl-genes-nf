import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "normalise_pairwise_tsv.py"
PAIRWISE_QC = Path(__file__).parents[1] / "bin" / "pairwise_qc.py"
BENCHMARK_PRODUCTS = Path(__file__).parents[1] / "bin" / "build_benchmark_products.py"
EXTRACT_BENCHMARK = Path(__file__).parents[1] / "bin" / "extract_pairwise_benchmark.py"
TESTDATA = Path(__file__).parents[1] / "testdata" / "pairwise"


class PairwiseContractTest(unittest.TestCase):
    def test_legacy_headers_are_normalised_to_source_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "legacy.tsv"
            output = root / "normalised.tsv"
            source.write_text(
                "ensembl_gene_id\tcat_gene_id\tens_to_cat_concordance_rate\n"
                "a1\tb1\t1.0\n",
                encoding="utf-8",
            )
            subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output)],
                check=True,
            )
            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines()[0],
                "source_a_gene_id\tsource_b_gene_id\tsource_a_to_source_b_concordance_rate",
            )

    def test_bed12_normalisation_is_structural_and_has_no_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "models.tsv"
            subprocess.run(
                [sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(TESTDATA / "source_a.bed12"),
                 "--source", "source_a", "--pair-id", "bed_pair", "--assembly", "fixture",
                 "--format", "bed12", "--output", str(output)],
                check=True,
            )
            lines = output.read_text(encoding="utf-8").splitlines()
            header = lines[0].split("\t")
            row = dict(zip(header, lines[1].split("\t")))
            self.assertIn("transcript_id", header)
            self.assertIn("gene_assignment_type", header)
            self.assertIn("locus_id", header)
            self.assertEqual(row["source_format"], "bed12")
            self.assertEqual(row["model_status"], "valid")
            self.assertEqual(row["original_gene_id"], "")
            self.assertEqual(row["gene_assignment_type"], "unassigned")
            self.assertEqual(row["cds_coordinates"], "")

    def test_bed12_conversion_marks_synthetic_gene_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "converted.gtf"
            subprocess.run(
                [sys.executable, str(PAIRWISE_QC), "bed12gtf", "--input", str(TESTDATA / "source_a.bed12"), "--output", str(output)],
                check=True,
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn('gene_id "synthetic_gene:lr_tx_1"', text)

    def test_manifest_validation_rejects_unknown_format(self):
        result = subprocess.run(
            [sys.executable, str(PAIRWISE_QC), "validate", "--pair-id", "p", "--assembly", "a",
             "--source-a", "a", "--source-b", "b", "--source-a-path", str(TESTDATA / "source_a.gff3"),
             "--source-b-path", str(TESTDATA / "source_b.gff3"), "--source-a-format", "wat",
             "--source-b-format", "gff3", "--output", "/tmp/invalid-pair.tsv"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_manifest_validation_rejects_mismatched_assemblies(self):
        result = subprocess.run(
            [sys.executable, str(PAIRWISE_QC), "validate", "--pair-id", "p", "--assembly", "a",
             "--source-a", "a", "--source-b", "b", "--source-a-path", str(TESTDATA / "source_a.gff3"),
             "--source-b-path", str(TESTDATA / "source_b.gff3"), "--source-a-format", "gff3",
             "--source-b-format", "gff3", "--source-a-assembly", "a", "--source-b-assembly", "b",
             "--output", "/tmp/invalid-assembly.tsv"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def _normalise(self, source, label, fmt, pair_id="p", assembly="fixture"):
        with tempfile.NamedTemporaryFile(suffix=".tsv") as output:
            subprocess.run(
                [sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(source),
                 "--source", label, "--pair-id", pair_id, "--assembly", assembly,
                 "--format", fmt, "--output", output.name], check=True,
            )
            return Path(output.name).read_text(encoding="utf-8")

    def test_gff3_source_gene_ids_are_preserved(self):
        text = self._normalise(TESTDATA / "source_a.gff3", "source_a", "gff3")
        row = dict(zip(text.splitlines()[0].split("\t"), text.splitlines()[1].split("\t")))
        self.assertEqual(row["original_gene_id"], "gene-a1")
        self.assertEqual(row["gene_assignment_type"], "biological_source_gene")

    def test_de_novo_loci_are_union_derived_and_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, b = root / "a.tsv", root / "b.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(TESTDATA / "source_a.bed12"), "--source", "source_a", "--pair-id", "p", "--assembly", "fixture", "--format", "bed12", "--output", str(a)], check=True)
            subprocess.run([sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(TESTDATA / "source_b.bed12"), "--source", "source_b", "--pair-id", "p", "--assembly", "fixture", "--format", "bed12", "--output", str(b)], check=True)
            outputs = []
            for suffix in ("one", "two"):
                out_a, out_b, summary = (root / f"a_{suffix}.tsv", root / f"b_{suffix}.tsv", root / f"summary_{suffix}.tsv")
                subprocess.run([sys.executable, str(PAIRWISE_QC), "assign-loci", "--source-a", str(a), "--source-b", str(b), "--pair-id", "p", "--assembly", "fixture", "--output-a", str(out_a), "--output-b", str(out_b), "--summary", str(summary)], check=True)
                outputs.append((out_a.read_text(), out_b.read_text()))
            self.assertEqual(outputs[0], outputs[1])
            self.assertIn("synthetic_locus", outputs[0][0])
            self.assertIn("locus:fixture:chr1:+:000001", outputs[0][0])
            self.assertIn("locus:fixture:chr1:+:000001", outputs[0][1])

    def test_locus_assignment_scales_with_many_independent_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = [
                "pair_id", "source", "assembly", "contig", "strand",
                "transcript_id", "original_transcript_id", "original_gene_id",
                "transcript_start", "transcript_end", "exon_count",
                "exon_coordinates", "intron_chain", "five_prime_end",
                "three_prime_end", "source_format", "gene_assignment_type",
                "locus_id", "locus_assignment_method", "model_status",
                "utr_coordinates", "cds_coordinates",
            ]

            def model(source, index):
                start = index * 1000 + 1
                end = start + 99
                transcript = f"{source}_{index}"
                return [
                    "p", source, "fixture", "chr1", "+", transcript, transcript,
                    "", str(start), str(end), "1", f"{start}-{end}", "",
                    str(start), str(end), "bed12", "unassigned", "", "", "valid",
                    "", "",
                ]
            for source in ("source_a", "source_b"):
                (root / f"{source}.tsv").write_text(
                    "\t".join(header) + "\n"
                    + "\n".join("\t".join(model(source, index)) for index in range(6000))
                    + "\n",
                    encoding="utf-8",
                )
            output_a, output_b, summary = (root / name for name in ("a.tsv", "b.tsv", "summary.tsv"))
            subprocess.run(
                [
                    sys.executable, str(PAIRWISE_QC), "assign-loci",
                    "--source-a", str(root / "source_a.tsv"),
                    "--source-b", str(root / "source_b.tsv"),
                    "--pair-id", "p", "--assembly", "fixture",
                    "--output-a", str(output_a), "--output-b", str(output_b),
                    "--summary", str(summary),
                ],
                check=True,
                timeout=20,
            )
            self.assertEqual(len(output_a.read_text(encoding="utf-8").splitlines()), 6001)
            self.assertEqual(len(output_b.read_text(encoding="utf-8").splitlines()), 6001)
            self.assertIn("synthetic_locus", summary.read_text(encoding="utf-8"))

    def test_reference_assignment_is_distinct_from_source_gene_assignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, b = root / "a.tsv", root / "b.tsv"
            for source, label, fmt, output in ((TESTDATA / "source_ref.gtf", "source_a", "gtf", a), (TESTDATA / "source_a.bed12", "source_b", "bed12", b)):
                subprocess.run([sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(source), "--source", label, "--pair-id", "p", "--assembly", "fixture", "--format", fmt, "--output", str(output)], check=True)
            out_a, out_b, summary = root / "ac.tsv", root / "bc.tsv", root / "summary.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "assign-loci", "--source-a", str(a), "--source-b", str(b), "--pair-id", "p", "--assembly", "fixture", "--reference", str(TESTDATA / "source_ref.gtf"), "--reference-format", "gtf", "--output-a", str(out_a), "--output-b", str(out_b), "--summary", str(summary)], check=True)
            row_a = dict(zip(out_a.read_text().splitlines()[0].split("\t"), out_a.read_text().splitlines()[1].split("\t")))
            row_b = dict(zip(out_b.read_text().splitlines()[0].split("\t"), out_b.read_text().splitlines()[1].split("\t")))
            self.assertEqual(row_a["gene_assignment_type"], "biological_source_gene")
            self.assertEqual(row_b["gene_assignment_type"], "reference_gene")
            self.assertEqual(row_b["locus_id"], "reference:gene1")
            summary_text = summary.read_text()
            self.assertRegex(summary_text, r"[0-9a-f]{64}")

    def test_synthetic_annotation_has_required_labels_and_no_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            models = root / "models.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(TESTDATA / "source_a.bed12"), "--source", "source_a", "--pair-id", "p", "--assembly", "fixture", "--format", "bed12", "--output", str(models)], check=True)
            assigned = root / "assigned.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "assign-loci", "--source-a", str(models), "--source-b", str(models), "--pair-id", "p", "--assembly", "fixture", "--output-a", str(assigned), "--output-b", str(root / "other.tsv"), "--summary", str(root / "summary.tsv")], check=True)
            annotation = root / "annotation.gff3"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "render-synthetic", "--input", str(assigned), "--source", "source_a", "--output", str(annotation)], check=True)
            text = annotation.read_text()
            self.assertIn("synthetic_gene=true", text)
            self.assertIn("source_transcript_id=lr_tx_1", text)
            self.assertNotIn("\tCDS\t", text)

    def test_metric_status_uses_explicit_capability_states(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "status.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "status", "--pair-id", "p", "--output", str(output), "--entry", "gene_presence", "not_applicable", "BED12", "--entry", "protein_validation", "missing_input", "no protein evidence"], check=True)
            self.assertIn("not_applicable", output.read_text())
            self.assertIn("missing_input", output.read_text())

    def test_metric_status_outputs_are_unique_for_run_level_aggregation(self):
        module = (Path(__file__).parents[1] / "modules" / "pairwise_annotation" / "metric_status.nf").read_text(encoding="utf-8")
        self.assertIn('path("${meta.pair_id}.metric_status.tsv"), emit: status', module)
        self.assertIn('> \'${meta.pair_id}.metric_status.tsv\'', module)
        self.assertNotIn("path('metric_status.tsv'), emit: status", module)

    def test_transcript_feature_comparison_separates_terminal_difference_from_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_a = root / "a.bed12"
            source_b = root / "b.bed12"
            source_a.write_text("chr1\t100\t220\ta_tx\t0\t+\t100\t220\t0\t2\t40,50\t0,70\n", encoding="utf-8")
            source_b.write_text("chr1\t100\t230\tb_tx\t0\t+\t100\t230\t0\t2\t40,60\t0,70\n", encoding="utf-8")
            a_models, b_models = root / "a.tsv", root / "b.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "normalise-pair", "--source-a", str(source_a), "--source-b", str(source_b), "--format-a", "bed12", "--format-b", "bed12", "--pair-id", "p", "--assembly", "fixture", "--output-a", str(a_models), "--output-b", str(b_models), "--audit", str(root / "audit.tsv")], check=True)
            output = root / "models.tsv"
            subprocess.run([sys.executable, str(PAIRWISE_QC), "compare-features", "--source-a", str(a_models), "--source-b", str(b_models), "--pair-id", "p", "--model-presence", str(output), "--feature-presence", str(root / "features.tsv"), "--terminal-exons", str(root / "terminal.tsv"), "--summary", str(root / "summary.tsv")], check=True)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["classification"] for row in rows}, {"same_intron_chain_3prime_terminal_difference"})
            self.assertEqual({row["same_intron_chain_models"] for row in rows}, {"1"})
            self.assertNotIn("locus_id", rows[0])

    def test_transcript_modules_use_one_command_without_printf_or_internal_versions(self):
        module_dir = Path(__file__).parents[1] / "modules" / "pairwise_annotation"
        for name in ("stage_transcript_inputs.nf", "normalise_transcript_pair.nf", "compare_transcript_features.nf"):
            text = (module_dir / name).read_text(encoding="utf-8")
            self.assertNotIn("printf", text)
            self.assertEqual(text.count("pairwise_qc.py"), 2)
            self.assertNotIn("versions.yml", text)

    def test_transcript_mode_is_explicitly_wired(self):
        main = (Path(__file__).parents[1] / "main.nf").read_text(encoding="utf-8")
        schema = (Path(__file__).parents[1] / "nextflow_schema.json").read_text(encoding="utf-8")
        subworkflow = (Path(__file__).parents[1] / "subworkflows" / "pairwise_annotation" / "pairwise_transcript_comparison.nf").read_text(encoding="utf-8")
        self.assertIn("pairwise_transcripts", main)
        self.assertIn("pairwise_transcripts", schema)
        self.assertNotIn("ASSIGN_TRANSCRIPT_LOCI", subworkflow)

    def test_benchmark_products_preserve_provenance_and_pairwise_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            models = {}
            for source, label, fmt in (("a", "bambu", "bed12"), ("b", "flair", "bed12")):
                models[source] = root / f"{source}.tsv"
                subprocess.run(
                    [sys.executable, str(PAIRWISE_QC), "normalise", "--input",
                     str(TESTDATA / f"source_{source}.bed12"), "--source", label,
                     "--pair-id", "p", "--assembly", "fixture", "--format", fmt,
                     "--output", str(models[source])], check=True,
                )
            metadata = root / "metadata.json"
            metadata.write_text(json.dumps({
                "pair_id": "p", "assembly": "fixture", "species": "human", "sample": "HG002",
                "source_a": "bambu", "source_b": "flair", "tool_a": "Bambu", "tool_b": "FLAIR",
                "technology_a": "PacBio_IsoSeq", "technology_b": "PacBio_IsoSeq",
                "format_a": "bed12", "format_b": "bed12",
            }), encoding="utf-8")
            intron = root / "intron.tsv"
            intron.write_text(
                "pair_id\tlocus_id\tsource_a_transcript_id\tsource_b_locus_id\tsource_b_transcript_id\tsource_a_introns\tsource_b_introns\tshared_introns\tjaccard\tclassification\n"
                "p\tlocus:fixture:chr1:+:000001\ttx_a\tlocus:fixture:chr1:+:000001\ttx_b\t1\t1\t1\t1.0\texact\n",
                encoding="utf-8",
            )
            ends = root / "ends.tsv"
            ends.write_text(
                "pair_id\tlocus_id\tend\tsource_a_transcript_id\tsource_b_locus_id\tsource_b_transcript_id\tdistance_bp\tclassification\n"
                "p\tlocus:fixture:chr1:+:000001\tfive_prime_end\ttx_a\tlocus:fixture:chr1:+:000001\ttx_b\t0\twithin_threshold\n"
                "p\tlocus:fixture:chr1:+:000001\tthree_prime_end\ttx_a\tlocus:fixture:chr1:+:000001\ttx_b\t0\twithin_threshold\n",
                encoding="utf-8",
            )
            utr = root / "utr.tsv"
            utr.write_text(
                "pair_id\tlocus_id\tsource_a_transcript_id\tsource_b_transcript_id\tsource_a_utr_available\tsource_b_utr_available\texact_match\tclassification\n"
                "p\tlocus:fixture:chr1:+:000001\ttx_a\ttx_b\t0\t0\t0\tmissing\n",
                encoding="utf-8",
            )
            structural = root / "structural.tsv"
            structural.write_text(
                "pair_id\tn_intron_comparisons\tn_exact_intron_chains\tn_partial_intron_chains\tn_end_comparisons\tn_ends_within_threshold\tn_utr_comparisons\tn_exact_utr_matches\n"
                "p\t1\t1\t0\t1\t1\t1\t0\n",
                encoding="utf-8",
            )
            output_models, output_comparisons, output_whole_models, output_summary = (root / name for name in ("models.tsv", "comparisons.tsv", "whole_models.tsv", "summary.tsv"))
            subprocess.run(
                [sys.executable, str(BENCHMARK_PRODUCTS), "pair", "--metadata", str(metadata),
                 "--source-a", str(models["a"]), "--source-b", str(models["b"]),
                 "--intron", str(intron), "--ends", str(ends), "--utr", str(utr),
                 "--structural", str(structural), "--output-models", str(output_models),
                 "--output-comparisons", str(output_comparisons), "--output-whole-models", str(output_whole_models),
                 "--output-summary", str(output_summary)],
                check=True,
            )
            self.assertIn("technology", output_models.read_text(encoding="utf-8").splitlines()[0])
            self.assertIn("Bambu", output_models.read_text(encoding="utf-8"))
            comparison_text = output_comparisons.read_text(encoding="utf-8")
            self.assertIn("intron_chain", comparison_text)
            self.assertIn("utr_structure", comparison_text)
            self.assertIn("whole_model_structure_only", output_whole_models.read_text(encoding="utf-8"))
            summary_text = output_summary.read_text(encoding="utf-8")
            self.assertIn("pairwise_concordance", summary_text)
            self.assertIn("none", summary_text)
            status = root / "status.tsv"
            status.write_text(
                "pair_id\tmetric\tstatus\treason\n"
                "p\tintron_chain\tenabled\tcanonical models\n"
                "p\tgene_presence\tnot_applicable\tBED12\n",
                encoding="utf-8",
            )
            subprocess.run(
                [sys.executable, str(BENCHMARK_PRODUCTS), "aggregate", "--models", str(output_models),
                 "--comparisons", str(output_comparisons), "--whole-models", str(output_whole_models),
                 "--summaries", str(output_summary), "--status", str(status)], cwd=root, check=True,
            )
            self.assertIn("Bambu", (root / "models_long.tsv").read_text(encoding="utf-8"))
            self.assertIn("whole_model_structure_only", (root / "whole_model_comparisons.tsv").read_text(encoding="utf-8"))
            self.assertIn("intron_chain_concordance", (root / "tool_technology_summary.tsv").read_text(encoding="utf-8"))
            report = (root / "benchmark_report.md").read_text(encoding="utf-8")
            self.assertIn("Long-read transcript reconstruction benchmark report", report)
            self.assertIn("Pairwise concordance is not recall", report)
            self.assertIn("not_applicable", (root / "metric_status_long.tsv").read_text(encoding="utf-8"))
            bundle = root / "bundle"
            subprocess.run(
                [sys.executable, str(EXTRACT_BENCHMARK), "--results-dir", str(root),
                 "--output-dir", str(bundle), "--pair-id", "p"], check=True,
            )
            self.assertIn("false", (bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue((bundle / "model_inventory_summary.tsv").is_file())
            self.assertTrue((bundle / "details" / "models.tsv.gz").is_file())


if __name__ == "__main__":
    unittest.main()
