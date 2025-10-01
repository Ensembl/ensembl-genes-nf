#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from collections import defaultdict
import hashlib


@dataclass
class Biogroup:
    id: str
    samples: List[str]
    shared_metadata: Dict
    sequencing_types: Set[str]
    study_sources: Set[str]


class BiogroupDetector:
    def __init__(self):
        # Semantic equivalences for biogroup detection
        self.treatment_equivalences = {
            "untreated": ["control", "mock", "vehicle", "dmso", "pbs"],
            "cycloheximide": ["chx", "cyclohexamide"],
            "harringtonine": ["hrr", "homoharringtonine"]
        }
        
        self.cell_type_equivalences = {
            "hek293": ["hek-293", "293", "hek293t"],
            "hela": ["he-la"],
            "cho": ["chinese hamster ovary"]
        }

    def normalize_condition(self, value: str) -> str:
        """Normalize a condition value for grouping."""
        if not value or value == "unknown":
            return "unknown"
            
        normalized = value.lower().strip()
        
        # Apply equivalence mappings
        for canonical, equivalents in self.treatment_equivalences.items():
            if normalized == canonical or normalized in equivalents:
                return canonical
                
        for canonical, equivalents in self.cell_type_equivalences.items():
            if normalized == canonical or normalized in equivalents:
                return canonical
        
        # Remove common units and normalize formatting
        normalized = normalized.replace("µm", "um").replace("μm", "um")
        normalized = normalized.replace(" ", "_")
        
        return normalized

    def create_biogroup_signature(self, sample_metadata: Dict) -> str:
        """Create a signature for biogroup clustering using standardized fields."""
        # Use standardized fields from evidence-based processing
        core_fields = ["cell_type", "treatment", "timepoint", "genotype", "tissue"]
        
        signature_parts = []
        for field in core_fields:
            standardized_field = f"standardized_{field}"
            
            if standardized_field in sample_metadata:
                # Use already standardized value
                value = sample_metadata[standardized_field]
                normalized_value = self.normalize_condition(str(value))
                signature_parts.append(f"{field}:{normalized_value}")
            else:
                # Fallback to old extraction method for backwards compatibility
                extracted = sample_metadata.get("extracted_metadata", {})
                if field in extracted:
                    field_value = extracted[field]
                    if isinstance(field_value, dict):
                        value = field_value.get("value", "unknown")
                    else:
                        value = field_value
                    
                    normalized_value = self.normalize_condition(str(value))
                    signature_parts.append(f"{field}:{normalized_value}")
                else:
                    signature_parts.append(f"{field}:unknown")
        
        signature = "|".join(signature_parts)
        
        # Create hash for shorter biogroup ID
        biogroup_hash = hashlib.md5(signature.encode()).hexdigest()[:8]
        return biogroup_hash

    def detect_biogroups(self, metadata: Dict) -> Dict:
        """Detect biogroups across all samples in the metadata."""
        biogroup_map = defaultdict(list)
        sample_to_biogroup = {}
        
        # Collect all samples across studies/series
        all_samples = []
        
        if "studies" in metadata:
            # SRA format
            for study in metadata["studies"]:
                study_id = study.get("bioproject_metadata", {}).get("bioproject_id", "unknown")
                for experiment in study.get("experiments", []):
                    for run in experiment.get("runs", []):
                        sample_info = {
                            "sample_id": run.get("run_id"),
                            "experiment_id": experiment.get("experiment_id"),
                            "study_id": study_id,
                            "sequencing_type": experiment.get("sequencing_type", "RNA-seq"),
                            "metadata": experiment
                        }
                        all_samples.append(sample_info)
                        
        elif "series" in metadata:
            # GEO format
            for series in metadata["series"]:
                series_id = series.get("series_metadata", {}).get("gse_id", "unknown")
                for sample in series.get("samples", []):
                    sample_info = {
                        "sample_id": sample.get("sample_id"),
                        "study_id": series_id,
                        "sequencing_type": sample.get("sequencing_type", "RNA-seq"),
                        "metadata": sample
                    }
                    all_samples.append(sample_info)
        
        # Group samples by biogroup signature
        for sample_info in all_samples:
            biogroup_id = self.create_biogroup_signature(sample_info["metadata"])
            biogroup_map[biogroup_id].append(sample_info)
            sample_to_biogroup[sample_info["sample_id"]] = biogroup_id
        
        # Create biogroup objects
        biogroups = []
        for biogroup_id, samples in biogroup_map.items():
            if len(samples) == 0:
                continue
                
            # Extract shared metadata
            shared_metadata = self._extract_shared_metadata(samples)
            
            # Collect sequencing types and study sources
            sequencing_types = set(sample["sequencing_type"] for sample in samples)
            study_sources = set(sample["study_id"] for sample in samples)
            
            biogroup = Biogroup(
                id=biogroup_id,
                samples=[sample["sample_id"] for sample in samples],
                shared_metadata=shared_metadata,
                sequencing_types=sequencing_types,
                study_sources=study_sources
            )
            
            biogroups.append(biogroup)
        
        # Sort biogroups by sample count (largest first)
        biogroups.sort(key=lambda bg: len(bg.samples), reverse=True)
        
        return {
            "biogroup_count": len(biogroups),
            "total_samples": len(all_samples),
            "cross_modal_biogroups": len([bg for bg in biogroups if len(bg.sequencing_types) > 1]),
            "cross_study_biogroups": len([bg for bg in biogroups if len(bg.study_sources) > 1]),
            "biogroups": [self._biogroup_to_dict(bg) for bg in biogroups],
            "sample_to_biogroup_map": sample_to_biogroup
        }

    def _extract_shared_metadata(self, samples: List[Dict]) -> Dict:
        """Extract metadata fields shared across all samples in a biogroup."""
        if not samples:
            return {}
        
        # Use standardized fields from evidence-based processing
        core_fields = ["cell_type", "treatment", "timepoint", "genotype", "tissue"]
        shared = {}
        
        for field in core_fields:
            standardized_field = f"standardized_{field}"
            field_values = []
            
            for sample in samples:
                sample_metadata = sample["metadata"]
                if standardized_field in sample_metadata:
                    field_values.append(sample_metadata[standardized_field])
                else:
                    # Fallback to extracted_metadata for backwards compatibility
                    extracted = sample_metadata.get("extracted_metadata", {})
                    if field in extracted:
                        field_value = extracted[field]
                        if isinstance(field_value, dict):
                            field_values.append(field_value.get("value", "unknown"))
                        else:
                            field_values.append(field_value)
                    else:
                        field_values.append("unknown")
            
            # If all values are the same, include in shared metadata
            unique_values = set(str(v) for v in field_values)
            if len(unique_values) == 1:
                shared[field] = field_values[0]
        
        return shared

    def _biogroup_to_dict(self, biogroup: Biogroup) -> Dict:
        """Convert biogroup object to dictionary for JSON serialization."""
        return {
            "biogroup_id": biogroup.id,
            "sample_count": len(biogroup.samples),
            "sample_ids": biogroup.samples,
            "shared_metadata": biogroup.shared_metadata,
            "sequencing_types": list(biogroup.sequencing_types),
            "study_sources": list(biogroup.study_sources),
            "is_cross_modal": len(biogroup.sequencing_types) > 1,
            "is_cross_study": len(biogroup.study_sources) > 1
        }

    def generate_biogroup_report(self, biogroup_data: Dict) -> Dict:
        """Generate summary report of biogroup detection results."""
        biogroups = biogroup_data["biogroups"]
        
        # Size distribution
        size_distribution = defaultdict(int)
        for bg in biogroups:
            sample_count = bg["sample_count"]
            if sample_count == 1:
                size_distribution["singleton"] += 1
            elif sample_count <= 5:
                size_distribution["small (2-5)"] += 1
            elif sample_count <= 20:
                size_distribution["medium (6-20)"] += 1
            else:
                size_distribution["large (20+)"] += 1
        
        # Cross-modal analysis
        cross_modal_groups = [bg for bg in biogroups if bg["is_cross_modal"]]
        cross_study_groups = [bg for bg in biogroups if bg["is_cross_study"]]
        
        # Most common biogroup patterns
        pattern_counts = defaultdict(int)
        for bg in biogroups:
            shared = bg["shared_metadata"]
            pattern_key = "|".join([
                f"cell_type:{shared.get('cell_type', 'unknown')}",
                f"treatment:{shared.get('treatment', 'unknown')}"
            ])
            pattern_counts[pattern_key] += 1
        
        top_patterns = sorted(pattern_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        return {
            "total_biogroups": len(biogroups),
            "singleton_samples": size_distribution["singleton"],
            "size_distribution": dict(size_distribution),
            "cross_modal_biogroups": len(cross_modal_groups),
            "cross_study_biogroups": len(cross_study_groups),
            "top_biogroup_patterns": top_patterns,
            "largest_biogroup_size": max((bg["sample_count"] for bg in biogroups), default=0)
        }


def main():
    parser = argparse.ArgumentParser(
        description="Detect biogroups from processed sample metadata"
    )
    parser.add_argument("input", help="Input JSON file with processed sample metadata")
    parser.add_argument("-o", "--output", help="Output JSON file with biogroup information", required=True)
    parser.add_argument("--report", help="Generate summary report file")
    
    args = parser.parse_args()
    
    detector = BiogroupDetector()
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    biogroup_data = detector.detect_biogroups(metadata)
    
    with open(args.output, 'w') as f:
        json.dump(biogroup_data, f, indent=2)
    
    if args.report:
        report = detector.generate_biogroup_report(biogroup_data)
        with open(args.report, 'w') as f:
            json.dump(report, f, indent=2)
    
    print(f"Detected {biogroup_data['biogroup_count']} biogroups from {biogroup_data['total_samples']} samples")
    print(f"Cross-modal biogroups: {biogroup_data['cross_modal_biogroups']}")
    print(f"Cross-study biogroups: {biogroup_data['cross_study_biogroups']}")


if __name__ == "__main__":
    main()