import csv
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import resolve_bam_reference
import validate_reference
import validate_star_junctions
import write_track_record


class BinHelperTests(unittest.TestCase):
    def test_resolver_keeps_bam_names_when_only_names_differ(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            chrom_sizes = root / "reference.chrom.sizes"
            chrom_sizes.write_text("CM001\t100\nNC_1\t20\n")
            header = "@HD\tVN:1.6\n@SQ\tSN:1\tLN:100\n@SQ\tSN:MT\tLN:20\n"
            results = [
                type("Completed", (), {"stdout": header})(),
                type("Completed", (), {"stdout": "samtools 1.20\n"})(),
            ]
            argv = ["resolve_bam_reference.py", "--bam", str(root / "x.bam"),
                    "--chrom-sizes", str(chrom_sizes), "--gca", "GCA_000000001.1",
                    "--assembly-release", "test", "--entity-id", "run_1",
                    "--output-prefix", str(root / "out")]
            with patch.object(sys, "argv", argv), patch.object(resolve_bam_reference.subprocess,
                    "run", side_effect=results), patch("sys.stderr", new_callable=io.StringIO):
                self.assertEqual(resolve_bam_reference.main(), 0)
            self.assertEqual((root / "out.chrom.sizes").read_text(), "1\t100\nMT\t20\n")
            with (root / "out.reference_resolution.tsv").open() as handle:
                report = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(report[1][3:], ["bam_names_used_lengths_match", "2"])

    def test_resolver_rejects_extra_bam_sequence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            chrom_sizes = root / "reference.chrom.sizes"
            chrom_sizes.write_text("CM001\t100\n")
            header = "@SQ\tSN:1\tLN:100\n@SQ\tSN:MT\tLN:20\n"
            argv = ["resolve_bam_reference.py", "--bam", str(root / "x.bam"),
                    "--chrom-sizes", str(chrom_sizes), "--gca", "GCA_000000001.1",
                    "--assembly-release", "test", "--entity-id", "run_1",
                    "--output-prefix", str(root / "out")]
            with patch.object(sys, "argv", argv), patch.object(resolve_bam_reference.subprocess,
                    "run", return_value=type("Completed", (), {"stdout": header})()), \
                    patch("sys.stderr", new_callable=io.StringIO) as stderr:
                self.assertEqual(resolve_bam_reference.main(), 1)
            self.assertIn("sequence lengths do not match", stderr.getvalue())

    def test_validate_reference_writes_same_check_contract(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            chrom_sizes = root / "chrom.sizes"
            chrom_sizes.write_text("chr1\t100\n")
            argv = ["validate_reference.py", "--chrom-sizes", str(chrom_sizes),
                    "--assembly-release", "test-release"]
            previous = Path.cwd()
            try:
                os.chdir(root)
                with patch.object(sys, "argv", argv):
                    self.assertEqual(validate_reference.main(), 0)
                self.assertIn("test-release", (root / "reference_check.tsv").read_text())
            finally:
                os.chdir(previous)

    def test_validate_junctions_accepts_header_and_data(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sj = root / "sample.SJ.out.tab"
            sj.write_text("chrom\tstart\tend\tstrand\tmotif\tannotated\tuniqReads\tmultiReads\tmaxOverhang\n"
                          "chr1\t2\t4\t1\t1\t0\t3\t1\t10\n")
            output = root / "validation.tsv"
            argv = ["validate_star_junctions.py", "--junctions", str(sj),
                    "--entity-id", "sample", "--output", str(output)]
            with patch.object(sys, "argv", argv):
                self.assertEqual(validate_star_junctions.main(), 0)
            self.assertIn("sample\tsplice_junction\tcomplete", output.read_text())

    def test_write_track_record_writes_result_and_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            track, checksum = root / "input.bw", root / "input.sha256"
            track.write_bytes(b"track")
            checksum.write_text("abc123  input.bw\n")
            artifact = root / "run.forward.bw"
            result = root / "run.result.tsv"
            provenance = root / "run.provenance.tsv"
            exceptions = root / "run.exceptions.tsv"
            versions = root / "run.versions.yml"
            argv = ["write_track_record.py", "--track", str(track), "--checksum", str(checksum),
                    "--artifact", str(artifact), "--track-result", str(result),
                    "--provenance", str(provenance), "--exceptions", str(exceptions),
                    "--versions", str(versions), "--gca", "GCA_000000001.1",
                    "--assembly-release", "test", "--entity-id", "run", "--entity-type", "run",
                    "--track-type", "coverage_forward", "--destination", "/out/run.forward.bw",
                    "--source-path", "/input/run.bam", "--normalization", '{"bin_size":25}']
            with patch.object(sys, "argv", argv):
                self.assertEqual(write_track_record.main(), 0)
            self.assertEqual(artifact.read_bytes(), b"track")
            result_text = result.read_text()
            self.assertIn("abc123", result_text)
            self.assertEqual(provenance.read_text(), result_text)
            self.assertEqual(exceptions.read_text(), "")
            self.assertIn("record_writer: pipeline", versions.read_text())


if __name__ == "__main__":
    unittest.main()
