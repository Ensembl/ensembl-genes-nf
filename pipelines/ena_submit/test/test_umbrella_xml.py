import csv
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
import unittest


PIPELINE_DIR = Path(__file__).parents[1]
GENERATOR = PIPELINE_DIR / "bin" / "generate_umbrella_xml.py"


class UmbrellaXmlTests(unittest.TestCase):
    def test_add_generates_umbrella_project(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir)
            command = [
                sys.executable,
                str(GENERATOR),
                "--action", "add",
                "--alias", "ensembl_umbrella",
                "--title", "Ensembl umbrella",
                "--description", "Umbrella for Ensembl evidence projects",
                "--outdir", str(output),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            tree = ET.parse(output / "umbrella_project.xml")
            self.assertIsNotNone(tree.find("./PROJECT_SET/PROJECT/UMBRELLA_PROJECT"))
            self.assertIsNotNone(tree.find("./SUBMISSION/ACTIONS/ACTION/ADD"))

    def test_modify_includes_successful_child_project_only(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            accessions = temp / "project_accessions.tsv"
            with accessions.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["analysis_id", "submission_id", "accession", "status"])
                writer.writerow(["child", "ERA1", "PRJEB123456", "SUCCESS"])
                writer.writerow(["failed", "ERA2", "", "FAILED"])

            command = [
                sys.executable,
                str(GENERATOR),
                "--action", "modify",
                "--alias", "ensembl_umbrella",
                "--title", "Ensembl umbrella",
                "--description", "Umbrella for Ensembl evidence projects",
                "--project-accessions", str(accessions),
                "--outdir", str(temp),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            tree = ET.parse(temp / "umbrella_project.xml")
            self.assertIsNotNone(tree.find("./SUBMISSION/ACTIONS/ACTION/MODIFY"))
            children = tree.findall("./PROJECT_SET/PROJECT/RELATED_PROJECTS/RELATED_PROJECT/CHILD_PROJECT")
            self.assertEqual([child.attrib["accession"] for child in children], ["PRJEB123456"])

    def test_dropbox_modify_separates_project_and_submission(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            project_accessions = temp / "project_accessions.tsv"
            analysis_accessions = temp / "analysis_accessions.tsv"
            project_accessions.write_text(
                "analysis_id\tsubmission_id\taccession\tstatus\n"
                "child\tERA1\tPRJEB123456\tSUCCESS\n"
            )
            analysis_accessions.write_text(
                "analysis_id\tsubmission_id\taccession\tstatus\n"
                "child\tERA2\tERZ123456\tSUCCESS\n"
            )
            command = [
                sys.executable,
                str(GENERATOR),
                "--action", "modify",
                "--alias", "ensembl_umbrella",
                "--title", "Ensembl umbrella",
                "--description", "Umbrella for Ensembl evidence projects",
                "--umbrella-accession", "PRJEB999999",
                "--format", "dropbox",
                "--project-accessions", str(project_accessions),
                "--analysis-accessions", str(analysis_accessions),
                "--outdir", str(temp),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            project = ET.parse(temp / "umbrella_project.xml")
            submission = ET.parse(temp / "umbrella_submission.xml")
            self.assertEqual(project.find("./PROJECT").attrib["accession"], "PRJEB999999")
            self.assertIsNotNone(project.find(".//CHILD_PROJECT[@accession='PRJEB123456']"))
            self.assertIsNotNone(submission.find("./ACTIONS/ACTION/MODIFY"))


if __name__ == "__main__":
    unittest.main()
