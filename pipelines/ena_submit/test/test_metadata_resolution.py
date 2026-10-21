import csv
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

from pipelines.ena_submit.bin import resolve_ena_metadata
from pipelines.ena_submit.bin import generate_analysis_xml
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
            "ERR1": {"run_accession": "ERR1", "sample_accession": "SAMEA1", "experiment_accession": "ERX1", "scientific_name": "Danio rerio", "tax_id": "7955"},
            "ERR2": {"run_accession": "ERR2", "sample_accession": "SAMEA2", "experiment_accession": "ERX2", "scientific_name": "Danio rerio", "tax_id": "7955"},
        }
        with patch.object(resolve_ena_metadata, "fetch_run", side_effect=lambda run, api: records[run]):
            resolved = resolve_ena_metadata.resolve_row(row, "unused")
        self.assertEqual(resolved["sample_accessions"], "SAMEA1,SAMEA2")
        self.assertEqual(resolved["experiment_accessions"], "ERX1,ERX2")
        self.assertEqual(resolved["taxon_id"], "7955")

    def test_biosample_secondary_accession_is_preferred(self):
        row = {"run_accession": "ERR1"}
        record = {
            "run_accession": "ERR1", "sample_accession": "ERS1",
            "secondary_sample_accession": "SAMEA1", "scientific_name": "Danio rerio", "tax_id": "7955",
        }
        with patch.object(resolve_ena_metadata, "fetch_run", return_value=record):
            resolved = resolve_ena_metadata.resolve_row(row, "unused")
        self.assertEqual(resolved["sample_accessions"], "SAMEA1")
        self.assertEqual(resolved["biosample_accessions"], "SAMEA1")

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
                "run_accessions": "ERR1", "sample_accessions": "SAMEA1",
                "scientific_name": "Danio rerio", "taxon_id": "7955",
                "aggregation_scope": "single_run", "project_alias": "prj_danio",
            },
            {
                "run_accessions": "ERR2", "sample_accessions": "SAMEA2",
                "scientific_name": "Danio rerio", "taxon_id": "7955",
                "aggregation_scope": "single_run", "project_alias": "prj_danio",
            },
        ]
        summary = resolve_ena_metadata.summarise_rows(rows)
        self.assertEqual(summary["sample_accessions"], "SAMEA1,SAMEA2")
        self.assertEqual(summary["scientific_name"], "Danio rerio")
        self.assertEqual(summary["taxon_id"], "7955")

    def test_analysis_emits_samples_and_schema_order(self):
        row = {
            "analysis_alias": "aggregate",
            "study": "ERP123456",
            "title": "Aggregate RNA-seq alignment evidence",
            "description": "Aggregate alignment evidence for a reference assembly.",
            "assembly_accession": "GCA_000001405.28",
            "sample_accessions": "SAMEA1,SAMN2",
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
        self.assertEqual([x.attrib["accession"] for x in analysis.findall("SAMPLE_REF")], ["SAMEA1", "SAMN2"])
        self.assertEqual(analysis.find("ANALYSIS_ATTRIBUTES/ANALYSIS_ATTRIBUTE/TAG").text, "species")

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
