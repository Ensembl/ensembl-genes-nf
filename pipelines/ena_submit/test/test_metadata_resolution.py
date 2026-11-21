import csv
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

from pipelines.ena_submit.bin import resolve_ena_metadata
from pipelines.ena_submit.bin import generate_analysis_xml
from pipelines.ena_submit.bin import generate_analysis_modify_xml
from pipelines.ena_submit.bin import generate_project_xml
from pipelines.ena_submit.bin import generate_umbrella_xml


class MetadataResolutionTests(unittest.TestCase):
    def test_aggregate_preserves_all_source_relationships(self):
        row = {
            "run_accessions": "ERR1,ERR2",
            "species": "Danio_rerio",
            "taxon_id": "7955",
        }
        records = {
            "ERR1": {"run_accession": "ERR1", "sample_accession": "ERS1", "secondary_sample_accession": "SAMEA1", "experiment_accession": "ERX1", "scientific_name": "Danio rerio", "tax_id": "7955"},
            "ERR2": {"run_accession": "ERR2", "sample_accession": "ERS2", "secondary_sample_accession": "SAMEA2", "experiment_accession": "ERX2", "scientific_name": "Danio rerio", "tax_id": "7955"},
        }
        with patch.object(resolve_ena_metadata, "fetch_run", side_effect=lambda run, api: records[run]):
            resolved = resolve_ena_metadata.resolve_row(row, "unused")
        self.assertEqual(resolved["sample_accessions"], "ERS1,ERS2")
        self.assertEqual(resolved["biosample_accessions"], "SAMEA1,SAMEA2")
        self.assertEqual(resolved["experiment_accessions"], "ERX1,ERX2")
        self.assertEqual(resolved["taxon_id"], "7955")

    def test_ena_sample_accession_is_used_and_biosample_is_retained(self):
        row = {"run_accession": "ERR1"}
        record = {
            "run_accession": "ERR1", "sample_accession": "ERS1",
            "secondary_sample_accession": "SAMEA1", "scientific_name": "Danio rerio", "tax_id": "7955",
        }
        with patch.object(resolve_ena_metadata, "fetch_run", return_value=record):
            resolved = resolve_ena_metadata.resolve_row(row, "unused")
        self.assertEqual(resolved["sample_accessions"], "ERS1")
        self.assertEqual(resolved["biosample_accessions"], "SAMEA1")

    def test_manifest_sample_mismatch_is_rejected(self):
        row = {"run_accession": "ERR1", "sample_accession": "SAMEA999"}
        record = {
            "run_accession": "ERR1", "sample_accession": "SAMEA1",
            "scientific_name": "Danio rerio", "tax_id": "7955",
        }
        with patch.object(resolve_ena_metadata, "fetch_run", return_value=record):
            with self.assertRaisesRegex(RuntimeError, "sample accession"):
                resolve_ena_metadata.resolve_row(row, "unused")

    def test_authoritative_metadata_replaces_manifest_values_and_is_retained(self):
        row = {
            "run_accession": "ERR1", "sample_accession": "ERS1",
            "experiment_accession": "ERX1", "platform": "ILLUMINA",
        }
        record = {
            "run_accession": "ERR1", "sample_accession": "ERS1",
            "secondary_sample_accession": "SAMEA1",
            "experiment_accession": "ERX1", "study_accession": "ERP1",
            "secondary_study_accession": "PRJEB1",
            "scientific_name": "Danio rerio", "tax_id": "7955",
            "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
            "instrument_platform": "ILLUMINA", "instrument_model": "Illumina NovaSeq 6000",
        }
        with patch.object(resolve_ena_metadata, "fetch_run", return_value=record):
            resolved = resolve_ena_metadata.resolve_row(row, "unused")
        self.assertEqual(resolved["sample_accession"], "ERS1")
        self.assertEqual(resolved["experiment_accession"], "ERX1")
        self.assertEqual(resolved["source_study_accessions"], "ERP1")
        self.assertEqual(resolved["secondary_study_accessions"], "PRJEB1")
        self.assertEqual(resolved["library_strategies"], "RNA-Seq")
        self.assertEqual(resolved["instrument_models"], "Illumina NovaSeq 6000")
        self.assertEqual(resolved["manifest_sample_accession"], "ERS1")

    def test_taxonomy_conflict_is_rejected(self):
        row = {"run_accessions": "ERR1,ERR2"}
        records = [
            {"run_accession": "ERR1", "scientific_name": "Danio rerio", "tax_id": "7955"},
            {"run_accession": "ERR2", "scientific_name": "Mus musculus", "tax_id": "10090"},
        ]
        with patch.object(resolve_ena_metadata, "fetch_run", side_effect=records):
            with self.assertRaisesRegex(RuntimeError, "disagree on taxonomy"):
                resolve_ena_metadata.resolve_row(row, "unused")

    def test_summary_collapses_resolved_metadata_for_project_creation(self):
        rows = [
            {
                "run_accessions": "ERR1", "sample_accessions": "ERS1",
                "scientific_name": "Danio rerio", "taxon_id": "7955",
                "aggregation_scope": "single_run", "project_alias": "prj_danio",
            },
            {
                "run_accessions": "ERR2", "sample_accessions": "ERS2",
                "scientific_name": "Danio rerio", "taxon_id": "7955",
                "aggregation_scope": "single_run", "project_alias": "prj_danio",
            },
        ]
        summary = resolve_ena_metadata.summarise_rows(rows)
        self.assertEqual(summary["sample_accessions"], "ERS1,ERS2")
        self.assertEqual(summary["scientific_name"], "Danio rerio")
        self.assertEqual(summary["taxon_id"], "7955")

    def test_analysis_emits_samples_and_schema_order(self):
        row = {
            "analysis_alias": "aggregate",
            "study": "ERP123456",
            "title": "Aggregate RNA-seq alignment evidence",
            "description": "Aggregate alignment evidence for a reference assembly.",
            "assembly_accession": "GCA_000001405.28",
            "sample_accessions": "ERS1,ERS2",
            "biosample_accessions": "SAMEA1,SAMN2",
            "experiment_accessions": "ERX1,ERX2",
            "run_accessions": "ERR1,ERR2",
            "species": "Homo_sapiens",
            "taxon_id": "9606",
        }
        files = [{"remote_path": "aggregate.bam", "file_type": "bam", "md5": "a" * 32}]
        root = generate_analysis_xml.build_analysis_xml(row, files)
        analysis = root.getroot().find("ANALYSIS")
        self.assertEqual([child.tag for child in analysis if child.tag.endswith("_REF")], [
            "STUDY_REF", "SAMPLE_REF", "SAMPLE_REF", "EXPERIMENT_REF", "EXPERIMENT_REF", "RUN_REF", "RUN_REF"
        ])
        self.assertEqual([x.attrib["accession"] for x in analysis.findall("SAMPLE_REF")], ["ERS1", "ERS2"])
        self.assertEqual(
            [ref.findtext("IDENTIFIERS/EXTERNAL_ID") for ref in analysis.findall("SAMPLE_REF")],
            ["SAMEA1", "SAMN2"],
        )
        attributes = {
            element.findtext("TAG"): element.findtext("VALUE")
            for element in analysis.findall("ANALYSIS_ATTRIBUTES/ANALYSIS_ATTRIBUTE")
        }
        self.assertEqual(attributes["species"], "Homo_sapiens")
        self.assertEqual(attributes["biosample_accessions"], "SAMEA1,SAMN2")

    def test_analysis_rejects_biosample_as_sample_ref(self):
        row = {
            "analysis_alias": "bad_sample",
            "study": "ERP123456",
            "assembly_accession": "GCA_000001405.28",
            "sample_accession": "SAMEA123456",
        }
        files = [{"remote_path": "alignment.cram", "file_type": "cram", "md5": "a" * 32}]
        with self.assertRaisesRegex(RuntimeError, "BioSample accession"):
            generate_analysis_xml.build_analysis_xml(row, files)

    def test_modify_resolves_ena_sample_and_preserves_file(self):
        root = ET.Element("ANALYSIS_SET")
        analysis = ET.SubElement(
            root,
            "ANALYSIS",
            accession="ERZ30663408",
            alias="rnaseq_alignment_evidence",
        )
        ET.SubElement(analysis, "STUDY_REF", accession="ERP209434")
        ET.SubElement(analysis, "RUN_REF", accession="ERR10476813")
        ET.SubElement(analysis, "SAMPLE_REF", accession="SAMEA111562625")
        analysis_type = ET.SubElement(analysis, "ANALYSIS_TYPE")
        ET.SubElement(analysis_type, "REFERENCE_ALIGNMENT")
        files = ET.SubElement(analysis, "FILES")
        ET.SubElement(
            files,
            "FILE",
            filename="analysis/ERZ306/ERZ30663408/alignment.bam",
            filetype="bam",
            checksum_method="MD5",
            checksum="a" * 32,
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "analysis.xml"
            ET.ElementTree(root).write(source, encoding="UTF-8", xml_declaration=True)
            record = {
                "run_accession": "ERR10476813",
                "sample_accession": "SAMEA111562625",
                "secondary_sample_accession": "ERS13672481",
            }
            with patch.object(generate_analysis_modify_xml, "fetch_run", return_value=record):
                tree, current = generate_analysis_modify_xml.read_current_analysis(source)
                samples, biosamples = generate_analysis_modify_xml.resolve_samples(
                    current, ["ERR10476813"], [], "unused"
                )
                generate_analysis_modify_xml.replace_sample_refs(current, samples, biosamples)
                generate_analysis_modify_xml.add_attribute(
                    current, "biosample_accessions", ",".join(biosamples)
                )
            modified = tree.getroot().find("ANALYSIS")
            self.assertEqual(
                [ref.attrib["accession"] for ref in modified.findall("SAMPLE_REF")],
                ["ERS13672481"],
            )
            self.assertEqual(
                modified.findtext("SAMPLE_REF/IDENTIFIERS/EXTERNAL_ID"),
                "SAMEA111562625",
            )
            self.assertEqual(
                modified.find("FILES/FILE").attrib["filename"],
                "analysis/ERZ306/ERZ30663408/alignment.bam",
            )
            self.assertEqual(
                modified.findtext("ANALYSIS_ATTRIBUTES/ANALYSIS_ATTRIBUTE/VALUE"),
                "SAMEA111562625",
            )
            self.assertIsNotNone(
                generate_analysis_modify_xml.build_submission_xml().getroot().find(
                    "ACTIONS/ACTION/MODIFY"
                )
            )
            self.assertIsNotNone(
                generate_analysis_modify_xml.build_submission_xml(validate=True).getroot().find(
                    "ACTIONS/ACTION/VALIDATE"
                )
            )

    def test_analysis_rejects_multiple_files_and_bad_md5(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            files = Path(temp_dir) / "files.tsv"
            files.write_text(
                "remote_path\tfile_type\tmd5\n"
                f"a.bam\tbam\t{'a' * 32}\n"
                f"b.bam\tbam\t{'b' * 32}\n"
            )
            with self.assertRaisesRegex(RuntimeError, "exactly one"):
                generate_analysis_xml.read_files_manifest(files)

            files.write_text("remote_path\tfile_type\tmd5\na.bam\tbam\tbad\n")
            with self.assertRaisesRegex(RuntimeError, "invalid MD5"):
                generate_analysis_xml.read_files_manifest(files)

    def test_project_uses_webin_v2_wrapper_and_structured_organism(self):
        root = generate_project_xml.build_webin_project(
            "prj_test", "prj_test", "Annotation evidence project for Danio rerio", 
            "A project containing reusable RNA-seq alignment evidence.",
            scientific_name="Danio rerio", taxon_id="7955"
        ).getroot()
        self.assertEqual(root.tag, "WEBIN")
        self.assertIsNotNone(root.find("SUBMISSION_SET/SUBMISSION/ACTIONS/ACTION/ADD"))
        organism = root.find("PROJECT_SET/PROJECT/SUBMISSION_PROJECT/ORGANISM")
        self.assertEqual(organism.findtext("TAXON_ID"), "7955")
        self.assertEqual(organism.findtext("SCIENTIFIC_NAME"), "Danio rerio")

    def test_umbrella_xml_adds_child_projects_with_modify(self):
        root = generate_umbrella_xml.build_xml(
            "PRJEB99999", ["PRJEB10001", "PRJEB10002"],
            "Ensembl alignment evidence umbrella study",
            "Umbrella study grouping Ensembl alignment evidence projects.",
        ).getroot()
        self.assertEqual(root.tag, "WEBIN")
        self.assertIsNotNone(root.find("SUBMISSION_SET/SUBMISSION/ACTIONS/ACTION/MODIFY"))
        children = root.findall("PROJECT_SET/PROJECT/RELATED_PROJECTS/RELATED_PROJECT/CHILD_PROJECT")
        self.assertEqual([child.attrib["accession"] for child in children], ["PRJEB10001", "PRJEB10002"])


if __name__ == "__main__":
    unittest.main()
