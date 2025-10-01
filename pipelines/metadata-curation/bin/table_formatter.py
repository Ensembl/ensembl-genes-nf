#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List, Optional, Any
import csv
from pathlib import Path


class TableFormatter:
    def __init__(self):
        # Define standard column schemas for different output types
        self.output_schemas = {
            "sample_level": [
                # Core identifiers
                "sample_id", "experiment_id", "study_id", 
                # Sequencing metadata
                "sequencing_type", "organism", "library_strategy", "library_source", "platform", "instrument",
                # Standardized fields (from evidence system)
                "standardized_tissue", "standardized_cell_type", "standardized_treatment", 
                "standardized_timepoint", "standardized_genotype", "standardized_replicate",
                # Run statistics
                "run_total_spots", "run_total_bases", "run_size", "run_published",
                # Study metadata
                "study_title", "study_description", "study_organism", "study_publication_pmid"
            ],
            "biogroup_level": [
                "biogroup_id", "sample_count", "sequencing_types", "study_sources",
                "cell_type", "treatment", "timepoint", "tissue", "organism",
                "is_cross_modal", "is_cross_study"
            ],
            "study_level": [
                "study_id", "title", "organism", "sample_count", "experiment_count",
                "sequencing_types", "publication_pmid", "submission_date"
            ]
        }

    def metadata_to_sample_table(self, metadata: Dict, include_raw: bool = False) -> List[Dict]:
        """Convert metadata to sample-level tabular format with complete SRA data plus standardized fields."""
        rows = []
        
        # Handle empty metadata
        if not metadata or (
            ("studies" in metadata and not metadata["studies"]) and 
            ("series" in metadata and not metadata["series"])
        ):
            print("Warning: No data to process, creating empty table", file=sys.stderr)
            return []
        
        if "studies" in metadata and metadata["studies"]:
            # SRA format - preserve ALL original data plus add standardized fields
            for study in metadata["studies"]:
                study_metadata = study.get("bioproject_metadata", {})
                study_id = study_metadata.get("bioproject_id")
                
                for experiment in study.get("experiments", []):
                    experiment_id = experiment.get("experiment_id")
                    
                    # Start with ALL experiment metadata (preserve everything)
                    base_row = {}
                    
                    # Add study-level metadata with prefix
                    for key, value in study_metadata.items():
                        base_row[f"study_{key}"] = value
                    
                    # Add experiment-level metadata
                    for key, value in experiment.items():
                        if key not in ["runs", "sample_attributes", "_evidence_standardization"]:
                            base_row[key] = value
                    
                    # Add ALL sample_attributes as individual columns
                    sample_attrs = experiment.get("sample_attributes", {})
                    for attr_name, attr_value in sample_attrs.items():
                        base_row[f"sra_{attr_name}"] = attr_value
                    
                    # Add standardized fields from evidence system
                    if "_evidence_standardization" in experiment:
                        evidence_data = experiment["_evidence_standardization"]
                        for field_name, field_info in evidence_data.get("standardized_fields", {}).items():
                            base_row[f"standardized_{field_name}"] = field_info.get("value", "unknown")
                            
                            # Include evidence metadata if requested
                            if include_raw and field_info.get("evidence_used"):
                                evidence = field_info["evidence_used"]
                                base_row[f"standardized_{field_name}_source"] = evidence.get("source")
                                base_row[f"standardized_{field_name}_method"] = evidence.get("method")
                                base_row[f"standardized_{field_name}_raw_value"] = evidence.get("raw_value")
                                base_row[f"standardized_{field_name}_evidence_count"] = field_info.get("all_evidence_count", 0)
                    
                    # Add legacy extracted metadata for backwards compatibility
                    extracted = experiment.get("extracted_metadata", {})
                    for field, value in extracted.items():
                        if isinstance(value, dict):
                            base_row[f"legacy_{field}"] = value.get("value")
                            if include_raw:
                                base_row[f"legacy_{field}_confidence"] = value.get("confidence")
                        else:
                            base_row[f"legacy_{field}"] = value
                    
                    # Create one row per run (this is the per-run granularity requested)
                    runs = experiment.get("runs", [])
                    if runs:
                        for run in runs:
                            row = base_row.copy()
                            # Add ALL run metadata
                            for run_key, run_value in run.items():
                                if run_key == "run_id":
                                    row["sample_id"] = run_value  # Primary identifier
                                row[f"run_{run_key}"] = run_value
                            rows.append(row)
                    else:
                        # No runs - use experiment as sample
                        row = base_row.copy()
                        row["sample_id"] = experiment.get("sample_id") or experiment_id
                        rows.append(row)
                        
        elif "series" in metadata and metadata["series"]:
            # GEO format
            for series in metadata["series"]:
                series_metadata = series.get("series_metadata", {})
                study_id = series_metadata.get("gse_id")
                
                for sample in series.get("samples", []):
                    row = {
                        "sample_id": sample.get("sample_id"),
                        "study_id": study_id,
                        "sequencing_type": sample.get("sequencing_type", "RNA-seq"),
                        "organism": sample.get("organism") or series_metadata.get("organism"),
                        "title": sample.get("title"),
                    }
                    
                    # Add extracted metadata
                    extracted = sample.get("extracted_metadata", {})
                    for field, value in extracted.items():
                        if isinstance(value, dict):
                            row[field] = value.get("value")
                            if include_raw:
                                row[f"{field}_confidence"] = value.get("confidence")
                        else:
                            row[field] = value
                    
                    rows.append(row)
        
        return rows

    def biogroups_to_table(self, biogroup_data: Dict) -> List[Dict]:
        """Convert biogroup data to tabular format."""
        rows = []
        
        for biogroup in biogroup_data.get("biogroups", []):
            row = {
                "biogroup_id": biogroup["biogroup_id"],
                "sample_count": biogroup["sample_count"],
                "sequencing_types": ",".join(biogroup["sequencing_types"]),
                "study_sources": ",".join(biogroup["study_sources"]),
                "is_cross_modal": biogroup["is_cross_modal"],
                "is_cross_study": biogroup["is_cross_study"]
            }
            
            # Add shared metadata fields
            shared_metadata = biogroup["shared_metadata"]
            for field, value in shared_metadata.items():
                row[field] = value
            
            rows.append(row)
        
        return rows

    def write_csv(self, rows: List[Dict], output_file: str, schema: Optional[List[str]] = None):
        """Write rows to CSV file with optional column ordering."""
        if not rows:
            print("No data to write, creating empty file with headers", file=sys.stderr)
            # Create empty file with just headers
            if schema:
                columns = schema
            else:
                columns = ["sample_id"]  # Minimal header
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=columns)
                writer.writeheader()
            return
        
        # Determine columns
        if schema:
            # Use provided schema, add any extra columns
            all_columns = set()
            for row in rows:
                all_columns.update(row.keys())
            
            columns = []
            # Add schema columns in order
            for col in schema:
                if col in all_columns:
                    columns.append(col)
            
            # Add remaining columns
            for col in sorted(all_columns):
                if col not in columns:
                    columns.append(col)
        else:
            # Use all columns found in data
            all_columns = set()
            for row in rows:
                all_columns.update(row.keys())
            columns = sorted(all_columns)
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    def write_tsv(self, rows: List[Dict], output_file: str, schema: Optional[List[str]] = None):
        """Write rows to TSV file."""
        if not rows:
            print("No data to write, creating empty file with headers", file=sys.stderr)
            # Create empty file with just headers
            if schema:
                columns = schema
            else:
                columns = ["sample_id"]  # Minimal header
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(columns)
            return
        
        # Similar to CSV but with tab delimiter
        if schema:
            all_columns = set()
            for row in rows:
                all_columns.update(row.keys())
            
            columns = []
            for col in schema:
                if col in all_columns:
                    columns.append(col)
            
            for col in sorted(all_columns):
                if col not in columns:
                    columns.append(col)
        else:
            all_columns = set()
            for row in rows:
                all_columns.update(row.keys())
            columns = sorted(all_columns)
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def format_metadata(self, metadata: Dict, biogroup_data: Optional[Dict] = None, 
                       output_format: str = "csv", output_level: str = "sample",
                       include_raw: bool = False) -> str:
        """Format metadata into tabular output."""
        
        if output_level == "sample":
            rows = self.metadata_to_sample_table(metadata, include_raw)
            schema = self.output_schemas["sample_level"]
        elif output_level == "biogroup" and biogroup_data:
            rows = self.biogroups_to_table(biogroup_data)
            schema = self.output_schemas["biogroup_level"]
        elif output_level == "study":
            rows = self.metadata_to_study_table(metadata)
            schema = self.output_schemas["study_level"]
        else:
            raise ValueError(f"Invalid output_level: {output_level}")
        
        return rows, schema

    def metadata_to_study_table(self, metadata: Dict) -> List[Dict]:
        """Convert metadata to study-level summary table."""
        rows = []
        
        if "studies" in metadata:
            # SRA format
            for study in metadata["studies"]:
                study_metadata = study.get("bioproject_metadata", {})
                experiments = study.get("experiments", [])
                
                sequencing_types = set()
                total_runs = 0
                
                for experiment in experiments:
                    seq_type = experiment.get("sequencing_type", "RNA-seq")
                    sequencing_types.add(seq_type)
                    total_runs += len(experiment.get("runs", []))
                
                row = {
                    "study_id": study_metadata.get("bioproject_id"),
                    "title": study_metadata.get("title"),
                    "organism": study_metadata.get("organism"),
                    "experiment_count": len(experiments),
                    "total_runs": total_runs,
                    "sequencing_types": ",".join(sorted(sequencing_types)),
                    "publication_pmid": study_metadata.get("publication_pmid"),
                    "submission_date": study_metadata.get("submission_date"),
                    "release_date": study_metadata.get("release_date")
                }
                
                rows.append(row)
                
        elif "series" in metadata:
            # GEO format
            for series in metadata["series"]:
                series_metadata = series.get("series_metadata", {})
                samples = series.get("samples", [])
                
                row = {
                    "study_id": series_metadata.get("gse_id"),
                    "title": series_metadata.get("title"),
                    "organism": series_metadata.get("organism"),
                    "sample_count": len(samples),
                    "platform_count": len(series_metadata.get("platform_ids", [])),
                    "submission_date": series_metadata.get("submission_date"),
                    "last_update_date": series_metadata.get("last_update_date")
                }
                
                if series_metadata.get("pubmed_ids"):
                    row["publication_pmid"] = ",".join(series_metadata["pubmed_ids"])
                
                rows.append(row)
        
        return rows


def main():
    parser = argparse.ArgumentParser(
        description="Format processed metadata into tabular outputs"
    )
    parser.add_argument("input", help="Input JSON file with processed metadata")
    parser.add_argument("-o", "--output", help="Output file path", required=True)
    parser.add_argument("--biogroups", help="Biogroup JSON file for biogroup-level output")
    parser.add_argument("--format", choices=["csv", "tsv", "json"], default="csv",
                        help="Output format")
    parser.add_argument("--level", choices=["sample", "biogroup", "study"], default="sample",
                        help="Output granularity level")
    parser.add_argument("--include-raw", action="store_true",
                        help="Include raw metadata and confidence scores")
    
    args = parser.parse_args()
    
    formatter = TableFormatter()
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    biogroup_data = None
    if args.biogroups:
        with open(args.biogroups) as f:
            biogroup_data = json.load(f)
    
    # Format the data
    rows, schema = formatter.format_metadata(
        metadata, 
        biogroup_data, 
        args.format, 
        args.level, 
        args.include_raw
    )
    
    # Write output
    if args.format == "csv":
        formatter.write_csv(rows, args.output, schema)
    elif args.format == "tsv":
        formatter.write_tsv(rows, args.output, schema)
    elif args.format == "json":
        with open(args.output, 'w') as f:
            json.dump(rows, f, indent=2)
    
    print(f"Formatted {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()