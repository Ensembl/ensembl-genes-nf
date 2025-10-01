#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List, Optional, Set
from dataclasses import dataclass
from collections import defaultdict, Counter


@dataclass
class ValidationResult:
    field: str
    issue_type: str
    severity: str  # "error", "warning", "info"
    message: str
    sample_id: Optional[str] = None
    study_id: Optional[str] = None


class QualityValidator:
    def __init__(self):
        # Minimum completeness thresholds per sequencing type
        self.completeness_thresholds = {
            "RNA-seq": {
                "required_fields": ["cell_type", "treatment"],
                "recommended_fields": ["timepoint", "replicate", "tissue"]
            },
            "Ribo-seq": {
                "required_fields": ["cell_type", "treatment", "translation_inhibitor"],
                "recommended_fields": ["fractionation_method", "timepoint", "replicate"]
            },
            "CAGE": {
                "required_fields": ["cell_type", "cap_selection_method"],
                "recommended_fields": ["treatment", "timepoint"]
            },
            "long-read": {
                "required_fields": ["platform", "library_prep"],
                "recommended_fields": ["read_length_dist", "size_selection"]
            }
        }
        
        # Known problematic values
        self.invalid_values = {
            "not applicable", "n/a", "na", "null", "none", "missing", 
            "not specified", "not provided", "unknown", ""
        }

    def validate_metadata_completeness(self, metadata: Dict) -> List[ValidationResult]:
        """Validate metadata completeness across all samples."""
        validation_results = []
        
        if "studies" in metadata:
            # SRA format
            for study in metadata["studies"]:
                study_id = study.get("bioproject_metadata", {}).get("bioproject_id")
                for experiment in study.get("experiments", []):
                    results = self._validate_experiment(experiment, study_id)
                    validation_results.extend(results)
                    
        elif "series" in metadata:
            # GEO format
            for series in metadata["series"]:
                series_id = series.get("series_metadata", {}).get("gse_id")
                for sample in series.get("samples", []):
                    results = self._validate_sample(sample, series_id)
                    validation_results.extend(results)
        
        return validation_results

    def _validate_experiment(self, experiment: Dict, study_id: str) -> List[ValidationResult]:
        """Validate a single experiment's metadata."""
        results = []
        
        experiment_id = experiment.get("experiment_id")
        sequencing_type = experiment.get("sequencing_type", "RNA-seq")
        extracted_metadata = experiment.get("extracted_metadata", {})
        
        # Check completeness
        schema = self.completeness_thresholds.get(sequencing_type, self.completeness_thresholds["RNA-seq"])
        
        # Required fields validation
        for field in schema["required_fields"]:
            if field not in extracted_metadata:
                results.append(ValidationResult(
                    field=field,
                    issue_type="missing_required",
                    severity="error",
                    message=f"Required field '{field}' missing for {sequencing_type}",
                    sample_id=experiment_id,
                    study_id=study_id
                ))
            else:
                value = self._extract_value(extracted_metadata[field])
                if self._is_invalid_value(value):
                    results.append(ValidationResult(
                        field=field,
                        issue_type="invalid_value",
                        severity="error",
                        message=f"Required field '{field}' has invalid value: '{value}'",
                        sample_id=experiment_id,
                        study_id=study_id
                    ))
        
        # Recommended fields validation
        for field in schema["recommended_fields"]:
            if field not in extracted_metadata:
                results.append(ValidationResult(
                    field=field,
                    issue_type="missing_recommended",
                    severity="warning",
                    message=f"Recommended field '{field}' missing for {sequencing_type}",
                    sample_id=experiment_id,
                    study_id=study_id
                ))
        
        # Sequencing-specific validation
        if sequencing_type == "Ribo-seq":
            results.extend(self._validate_ribo_seq(extracted_metadata, experiment_id, study_id))
        elif sequencing_type == "CAGE":
            results.extend(self._validate_cage(extracted_metadata, experiment_id, study_id))
        
        return results

    def _validate_sample(self, sample: Dict, study_id: str) -> List[ValidationResult]:
        """Validate a single sample's metadata (GEO format)."""
        results = []
        
        sample_id = sample.get("sample_id")
        sequencing_type = sample.get("sequencing_type", "RNA-seq")
        extracted_metadata = sample.get("extracted_metadata", {})
        
        schema = self.completeness_thresholds.get(sequencing_type, self.completeness_thresholds["RNA-seq"])
        
        # Similar validation as for experiments
        for field in schema["required_fields"]:
            if field not in extracted_metadata:
                results.append(ValidationResult(
                    field=field,
                    issue_type="missing_required",
                    severity="error",
                    message=f"Required field '{field}' missing for {sequencing_type}",
                    sample_id=sample_id,
                    study_id=study_id
                ))
        
        return results

    def _validate_ribo_seq(self, metadata: Dict, sample_id: str, study_id: str) -> List[ValidationResult]:
        """Validate Ribo-seq specific requirements."""
        results = []
        
        # Check translation inhibitor
        inhibitor = self._extract_value(metadata.get("translation_inhibitor"))
        known_inhibitors = ["cycloheximide", "harringtonine", "puromycin", "chloramphenicol"]
        
        if inhibitor and inhibitor.lower() not in [inh.lower() for inh in known_inhibitors]:
            results.append(ValidationResult(
                field="translation_inhibitor",
                issue_type="unknown_inhibitor",
                severity="warning",
                message=f"Unknown translation inhibitor: '{inhibitor}'. Common inhibitors: {', '.join(known_inhibitors)}",
                sample_id=sample_id,
                study_id=study_id
            ))
        
        return results

    def _validate_cage(self, metadata: Dict, sample_id: str, study_id: str) -> List[ValidationResult]:
        """Validate CAGE-specific requirements."""
        results = []
        
        # Check cap selection method
        cap_method = self._extract_value(metadata.get("cap_selection_method"))
        known_methods = ["cap trapper", "oligo capping", "cap-seq"]
        
        if cap_method and cap_method.lower() not in [method.lower() for method in known_methods]:
            results.append(ValidationResult(
                field="cap_selection_method",
                issue_type="unknown_method",
                severity="warning",
                message=f"Unknown cap selection method: '{cap_method}'. Common methods: {', '.join(known_methods)}",
                sample_id=sample_id,
                study_id=study_id
            ))
        
        return results

    def _extract_value(self, field_data) -> str:
        """Extract value from field data (handles both string and dict formats)."""
        if isinstance(field_data, dict):
            return str(field_data.get("value", ""))
        return str(field_data) if field_data is not None else ""

    def _is_invalid_value(self, value: str) -> bool:
        """Check if a value is considered invalid."""
        return value.lower().strip() in self.invalid_values

    def validate_ontology_mappings(self, metadata: Dict) -> List[ValidationResult]:
        """Validate ontology mapping quality."""
        results = []
        
        # Check mapping results if available
        if "_mapping_results" in metadata:
            mapping_results = metadata["_mapping_results"]
            
            unmapped_count = sum(1 for result in mapping_results.values() 
                               if result.get("confidence") == "unmapped")
            total_mappings = len(mapping_results)
            
            if total_mappings > 0:
                unmapped_rate = unmapped_count / total_mappings
                if unmapped_rate > 0.5:
                    results.append(ValidationResult(
                        field="ontology_mappings",
                        issue_type="low_mapping_rate",
                        severity="warning",
                        message=f"High unmapped rate: {unmapped_rate:.2%} ({unmapped_count}/{total_mappings})"
                    ))
        
        return results

    def validate_biogroups(self, biogroup_data: Dict) -> List[ValidationResult]:
        """Validate biogroup detection quality."""
        results = []
        
        biogroups = biogroup_data.get("biogroups", [])
        total_samples = biogroup_data.get("total_samples", 0)
        
        # Check for too many singleton biogroups
        singletons = len([bg for bg in biogroups if bg["sample_count"] == 1])
        if total_samples > 0:
            singleton_rate = singletons / total_samples
            if singleton_rate > 0.7:
                results.append(ValidationResult(
                    field="biogroups",
                    issue_type="high_singleton_rate",
                    severity="warning",
                    message=f"High singleton rate: {singleton_rate:.2%} ({singletons}/{total_samples}). May indicate poor metadata quality."
                ))
        
        # Check for very large biogroups (might indicate over-grouping)
        large_groups = [bg for bg in biogroups if bg["sample_count"] > 50]
        if large_groups:
            results.append(ValidationResult(
                field="biogroups",
                issue_type="oversized_biogroups",
                severity="info",
                message=f"Found {len(large_groups)} biogroups with >50 samples. Verify grouping logic."
            ))
        
        return results

    def generate_validation_report(self, metadata: Dict, biogroup_data: Optional[Dict] = None) -> Dict:
        """Generate comprehensive validation report."""
        all_results = []
        
        # Validate metadata completeness
        all_results.extend(self.validate_metadata_completeness(metadata))
        
        # Validate ontology mappings
        all_results.extend(self.validate_ontology_mappings(metadata))
        
        # Validate biogroups if provided
        if biogroup_data:
            all_results.extend(self.validate_biogroups(biogroup_data))
        
        # Summarize results
        error_count = len([r for r in all_results if r.severity == "error"])
        warning_count = len([r for r in all_results if r.severity == "warning"])
        info_count = len([r for r in all_results if r.severity == "info"])
        
        # Group issues by type
        issues_by_type = defaultdict(list)
        for result in all_results:
            issues_by_type[result.issue_type].append(result)
        
        # Count samples with issues
        samples_with_errors = len(set(r.sample_id for r in all_results if r.sample_id and r.severity == "error"))
        
        report = {
            "validation_summary": {
                "total_issues": len(all_results),
                "errors": error_count,
                "warnings": warning_count,
                "info": info_count,
                "samples_with_errors": samples_with_errors
            },
            "issues_by_type": {
                issue_type: [
                    {
                        "field": r.field,
                        "severity": r.severity,
                        "message": r.message,
                        "sample_id": r.sample_id,
                        "study_id": r.study_id
                    }
                    for r in results
                ]
                for issue_type, results in issues_by_type.items()
            },
            "field_completion_rates": self._calculate_field_completion_rates(metadata),
            "data_quality_score": self._calculate_quality_score(all_results, metadata)
        }
        
        return report

    def _calculate_field_completion_rates(self, metadata: Dict) -> Dict:
        """Calculate completion rates for each metadata field."""
        field_counts = defaultdict(int)
        total_samples = 0
        
        if "studies" in metadata:
            for study in metadata["studies"]:
                for experiment in study.get("experiments", []):
                    total_samples += 1
                    extracted = experiment.get("extracted_metadata", {})
                    for field, value in extracted.items():
                        field_value = self._extract_value(value)
                        if not self._is_invalid_value(field_value):
                            field_counts[field] += 1
                            
        elif "series" in metadata:
            for series in metadata["series"]:
                for sample in series.get("samples", []):
                    total_samples += 1
                    extracted = sample.get("extracted_metadata", {})
                    for field, value in extracted.items():
                        field_value = self._extract_value(value)
                        if not self._is_invalid_value(field_value):
                            field_counts[field] += 1
        
        completion_rates = {}
        for field, count in field_counts.items():
            completion_rates[field] = count / total_samples if total_samples > 0 else 0
            
        return completion_rates

    def _calculate_quality_score(self, validation_results: List[ValidationResult], metadata: Dict) -> float:
        """Calculate overall data quality score (0-1)."""
        if not validation_results:
            return 1.0
        
        # Count total samples
        total_samples = 0
        if "studies" in metadata:
            total_samples = sum(len(study.get("experiments", [])) for study in metadata["studies"])
        elif "series" in metadata:
            total_samples = sum(len(series.get("samples", [])) for series in metadata["series"])
        
        if total_samples == 0:
            return 0.0
        
        # Calculate penalty weights
        error_penalty = 1.0
        warning_penalty = 0.5
        info_penalty = 0.1
        
        total_penalty = 0
        for result in validation_results:
            if result.severity == "error":
                total_penalty += error_penalty
            elif result.severity == "warning":
                total_penalty += warning_penalty
            elif result.severity == "info":
                total_penalty += info_penalty
        
        # Normalize by total samples
        normalized_penalty = total_penalty / total_samples
        
        # Convert to score (1 = perfect, 0 = very poor)
        quality_score = max(0.0, 1.0 - normalized_penalty)
        
        return quality_score

    def _extract_value(self, field_data) -> str:
        """Extract value from field data."""
        if isinstance(field_data, dict):
            return str(field_data.get("value", ""))
        return str(field_data) if field_data is not None else ""

    def _is_invalid_value(self, value: str) -> bool:
        """Check if a value is considered invalid."""
        return value.lower().strip() in self.invalid_values

    def validate_biogroup_quality(self, biogroup_data: Dict) -> List[ValidationResult]:
        """Validate biogroup detection quality."""
        results = []
        
        biogroups = biogroup_data.get("biogroups", [])
        
        # Check for reasonable biogroup sizes
        for biogroup in biogroups:
            sample_count = biogroup["sample_count"]
            biogroup_id = biogroup["biogroup_id"]
            
            # Flag very small biogroups (might indicate over-splitting)
            if sample_count == 1:
                results.append(ValidationResult(
                    field="biogroup_size",
                    issue_type="singleton_biogroup",
                    severity="info",
                    message=f"Biogroup {biogroup_id} contains only 1 sample",
                ))
            
            # Flag very large biogroups (might indicate under-splitting)
            elif sample_count > 100:
                results.append(ValidationResult(
                    field="biogroup_size",
                    issue_type="oversized_biogroup", 
                    severity="warning",
                    message=f"Biogroup {biogroup_id} contains {sample_count} samples - verify grouping logic",
                ))
        
        return results

    def validate_cross_modal_consistency(self, biogroup_data: Dict, metadata: Dict) -> List[ValidationResult]:
        """Validate consistency across different sequencing modalities."""
        results = []
        
        cross_modal_biogroups = [bg for bg in biogroup_data.get("biogroups", []) 
                               if bg.get("is_cross_modal", False)]
        
        for biogroup in cross_modal_biogroups:
            # Check that cross-modal samples have consistent core metadata
            biogroup_id = biogroup["biogroup_id"]
            sequencing_types = biogroup["sequencing_types"]
            
            # Ensure essential fields are consistent
            shared_metadata = biogroup["shared_metadata"]
            essential_fields = ["cell_type", "treatment"]
            
            for field in essential_fields:
                if field not in shared_metadata or self._is_invalid_value(str(shared_metadata[field])):
                    results.append(ValidationResult(
                        field=field,
                        issue_type="cross_modal_inconsistency",
                        severity="error",
                        message=f"Cross-modal biogroup {biogroup_id} missing consistent {field} across {', '.join(sequencing_types)}",
                    ))
        
        return results


def main():
    parser = argparse.ArgumentParser(
        description="Validate metadata quality and completeness"
    )
    parser.add_argument("input", help="Input JSON file with processed metadata")
    parser.add_argument("-o", "--output", help="Output JSON file with validation report", required=True)
    parser.add_argument("--biogroups", help="Biogroup JSON file for cross-validation")
    parser.add_argument("--summary-only", action="store_true",
                        help="Output only summary statistics")
    
    args = parser.parse_args()
    
    validator = QualityValidator()
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    biogroup_data = None
    if args.biogroups:
        with open(args.biogroups) as f:
            biogroup_data = json.load(f)
    
    # Generate validation report
    validation_report = validator.generate_validation_report(metadata, biogroup_data)
    
    # Add biogroup-specific validations if available
    if biogroup_data:
        biogroup_validation = validator.validate_biogroup_quality(biogroup_data)
        cross_modal_validation = validator.validate_cross_modal_consistency(biogroup_data, metadata)
        
        all_biogroup_issues = biogroup_validation + cross_modal_validation
        validation_report["biogroup_validation"] = [
            {
                "field": r.field,
                "issue_type": r.issue_type,
                "severity": r.severity,
                "message": r.message,
                "sample_id": r.sample_id,
                "study_id": r.study_id
            }
            for r in all_biogroup_issues
        ]
    
    if args.summary_only:
        summary_report = {
            "validation_summary": validation_report["validation_summary"],
            "data_quality_score": validation_report["data_quality_score"],
            "field_completion_rates": validation_report["field_completion_rates"]
        }
        output_report = summary_report
    else:
        output_report = validation_report
    
    with open(args.output, 'w') as f:
        json.dump(output_report, f, indent=2)
    
    summary = validation_report["validation_summary"]
    print(f"Validation completed. Quality score: {validation_report['data_quality_score']:.2f}")
    print(f"Issues found: {summary['errors']} errors, {summary['warnings']} warnings, {summary['info']} info")


if __name__ == "__main__":
    main()