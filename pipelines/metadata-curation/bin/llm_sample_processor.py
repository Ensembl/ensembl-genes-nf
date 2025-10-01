#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List, Optional, Any
import re


class LLMSampleProcessor:
    def __init__(self, model_endpoint: Optional[str] = None, api_key: Optional[str] = None):
        self.model_endpoint = model_endpoint
        self.api_key = api_key
        
        # Schema definitions for different sequencing types
        self.sequencing_schemas = {
            "RNA-seq": {
                "required_fields": ["cell_type", "treatment", "timepoint", "replicate"],
                "optional_fields": ["tissue", "development_stage", "genotype", "strain"]
            },
            "Ribo-seq": {
                "required_fields": ["translation_inhibitor", "fractionation_method", "cell_type", "treatment"],
                "optional_fields": ["chase_time", "drug_concentration", "ribosome_fraction"]
            },
            "CAGE": {
                "required_fields": ["cap_selection_method", "cell_type", "treatment"],
                "optional_fields": ["tss_enrichment", "degradome_treatment"]
            },
            "long-read": {
                "required_fields": ["platform", "library_prep", "read_length_dist"],
                "optional_fields": ["size_selection", "polymerase", "flow_cell_type"]
            }
        }

    def detect_sequencing_type(self, experiment_metadata: Dict) -> str:
        """Detect sequencing type from experiment metadata."""
        library_strategy = (experiment_metadata.get("library_strategy") or "").lower()
        library_source = (experiment_metadata.get("library_source") or "").lower()
        platform = (experiment_metadata.get("platform") or "").lower()
        title = (experiment_metadata.get("experiment_title") or "").lower()
        
        # Ribo-seq detection
        ribo_keywords = ["ribosome profiling", "ribo-seq", "ribo seq", "ribosome footprint"]
        if any(keyword in title or keyword in library_strategy for keyword in ribo_keywords):
            return "Ribo-seq"
        
        # CAGE detection  
        cage_keywords = ["cage", "cap analysis", "tss", "transcription start"]
        if any(keyword in title or keyword in library_strategy for keyword in cage_keywords):
            return "CAGE"
        
        # Long-read detection
        long_read_platforms = ["pacbio", "nanopore", "oxford", "sequel", "minion"]
        if any(platform_name in platform for platform_name in long_read_platforms):
            return "long-read"
        
        # Default to RNA-seq for transcriptomic data
        if library_source in ["transcriptomic", "genomic"] or "rna" in library_strategy:
            return "RNA-seq"
            
        return "RNA-seq"  # Default fallback

    def extract_sample_fields(self, sample_metadata: Dict, study_context: Dict, 
                            sequencing_type: str) -> Dict:
        """Extract and standardize sample-level fields using study context."""
        schema = self.sequencing_schemas.get(sequencing_type, self.sequencing_schemas["RNA-seq"])
        extracted_fields = {}
        
        # Get raw sample attributes
        sample_attrs = sample_metadata.get("sample_attributes", {})
        characteristics = sample_metadata.get("characteristics", {})
        
        # Combine all possible sources of sample info
        all_attributes = {**sample_attrs, **characteristics}
        
        # Add title and description as searchable text
        searchable_text = " ".join([
            sample_metadata.get("sample_title") or "",
            sample_metadata.get("title") or "",
            sample_metadata.get("description") or "",
            sample_metadata.get("treatment_protocol") or "",
            sample_metadata.get("growth_protocol") or ""
        ]).lower()
        
        # Extract fields using study context and pattern matching
        extracted_fields.update(self._extract_cell_type(all_attributes, searchable_text, study_context))
        extracted_fields.update(self._extract_treatment(all_attributes, searchable_text, study_context))
        extracted_fields.update(self._extract_timepoint(all_attributes, searchable_text, study_context))
        extracted_fields.update(self._extract_replicate(all_attributes, searchable_text))
        
        # Sequencing type specific extractions
        if sequencing_type == "Ribo-seq":
            extracted_fields.update(self._extract_ribo_seq_fields(all_attributes, searchable_text))
        elif sequencing_type == "CAGE":
            extracted_fields.update(self._extract_cage_fields(all_attributes, searchable_text))
        elif sequencing_type == "long-read":
            extracted_fields.update(self._extract_long_read_fields(all_attributes, searchable_text))
        
        # Add confidence scores
        for field, value in extracted_fields.items():
            if isinstance(value, dict) and "value" in value:
                continue  # Already has confidence info
            extracted_fields[field] = {
                "value": value,
                "confidence": self._calculate_confidence(field, value, all_attributes),
                "extraction_method": "pattern_matching"
            }
        
        return extracted_fields

    def _extract_cell_type(self, attributes: Dict, text: str, context: Dict) -> Dict:
        """Extract cell type information."""
        # Look in attributes first
        cell_type_keys = ["cell type", "cell_type", "celltype", "cell line", "cell_line"]
        for key in cell_type_keys:
            if key in attributes:
                return {"cell_type": attributes[key]}
        
        # Pattern matching in text
        cell_patterns = [
            r'\b([A-Z]{2,}[-0-9]*)\s*cells?\b',  # HEK293, CHO, etc.
            r'\b(primary\s+[^,\s]+)\s*cells?\b',
            r'\b([a-z]+\s+fibroblasts?|hepatocytes?|neurons?)\b'
        ]
        
        for pattern in cell_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                return {"cell_type": match.group(1)}
        
        # Use study context expectations
        expected_cells = context.get("expected_cell_types", [])
        if expected_cells:
            return {"cell_type": expected_cells[0]}  # Use first expected
            
        return {"cell_type": "unknown"}

    def _extract_treatment(self, attributes: Dict, text: str, context: Dict) -> Dict:
        """Extract treatment information."""
        treatment_keys = ["treatment", "compound", "drug", "stimulus", "condition"]
        for key in treatment_keys:
            if key in attributes:
                return {"treatment": attributes[key]}
        
        # Pattern matching
        treatment_patterns = [
            r'treated with ([^,\.]+)',
            r'(\d+\s*(?:µM|uM|mM|nM))\s+([A-Za-z][A-Za-z0-9\-]+)',
            r'\b(control|untreated|mock|vehicle)\b'
        ]
        
        for pattern in treatment_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                return {"treatment": match.group(1) if len(match.groups()) == 1 else match.group(2)}
        
        return {"treatment": "unknown"}

    def _extract_timepoint(self, attributes: Dict, text: str, context: Dict) -> Dict:
        """Extract timepoint information."""
        time_keys = ["time", "timepoint", "time point", "harvest time", "collection time"]
        for key in time_keys:
            if key in attributes:
                return {"timepoint": attributes[key]}
        
        # Pattern matching for time expressions
        time_patterns = [
            r'\b(\d+(?:\.\d+)?)\s*(hours?|hrs?|h)\b',
            r'\b(\d+(?:\.\d+)?)\s*(days?|d)\b',
            r'\b(\d+(?:\.\d+)?)\s*(minutes?|min)\b'
        ]
        
        for pattern in time_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                return {"timepoint": f"{match.group(1)} {match.group(2)}"}
        
        return {"timepoint": "unknown"}

    def _extract_replicate(self, attributes: Dict, text: str) -> Dict:
        """Extract replicate information."""
        rep_keys = ["replicate", "rep", "biological replicate", "biol rep"]
        for key in rep_keys:
            if key in attributes:
                return {"replicate": attributes[key]}
        
        # Pattern matching for replicate numbers
        rep_patterns = [
            r'\b(?:rep|replicate)\s*(\d+)\b',
            r'\b(\d+)(?:st|nd|rd|th)?\s*(?:rep|replicate)\b',
            r'\b(rep[A-D]|[A-D]rep)\b'
        ]
        
        for pattern in rep_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                return {"replicate": match.group(1)}
        
        return {"replicate": "unknown"}

    def _extract_ribo_seq_fields(self, attributes: Dict, text: str) -> Dict:
        """Extract Ribo-seq specific fields."""
        fields = {}
        
        # Translation inhibitor
        inhibitor_keys = ["inhibitor", "translation inhibitor", "drug"]
        for key in inhibitor_keys:
            if key in attributes:
                fields["translation_inhibitor"] = attributes[key]
                break
        
        if "translation_inhibitor" not in fields:
            inhibitors = ["cycloheximide", "harringtonine", "puromycin", "chloramphenicol"]
            for inhibitor in inhibitors:
                if inhibitor in text:
                    fields["translation_inhibitor"] = inhibitor
                    break
            else:
                fields["translation_inhibitor"] = "unknown"
        
        # Fractionation method
        if "monosome" in text or "polysome" in text:
            fields["fractionation_method"] = "sucrose gradient"
        else:
            fields["fractionation_method"] = "unknown"
            
        return fields

    def _extract_cage_fields(self, attributes: Dict, text: str) -> Dict:
        """Extract CAGE-specific fields."""
        fields = {}
        
        # Cap selection method
        if "cap trapper" in text or "cap-trapper" in text:
            fields["cap_selection_method"] = "cap trapper"
        elif "oligo capping" in text or "oligo-capping" in text:
            fields["cap_selection_method"] = "oligo capping"
        else:
            fields["cap_selection_method"] = "unknown"
            
        return fields

    def _extract_long_read_fields(self, attributes: Dict, text: str) -> Dict:
        """Extract long-read sequencing specific fields."""
        fields = {}
        
        # Platform detection
        if "pacbio" in text or "sequel" in text:
            fields["platform"] = "PacBio"
        elif "nanopore" in text or "minion" in text or "gridion" in text:
            fields["platform"] = "Oxford Nanopore"
        else:
            fields["platform"] = "unknown"
            
        return fields

    def _calculate_confidence(self, field: str, value: str, attributes: Dict) -> str:
        """Calculate confidence score for extracted field."""
        if value == "unknown":
            return "low"
        
        # High confidence if found in structured attributes
        if any(field.lower() in key.lower() for key in attributes.keys()):
            return "high"
        
        # Medium confidence for pattern matching
        return "medium"

    def process_samples_with_context(self, metadata: Dict) -> Dict:
        """Process all samples with study context."""
        processed_metadata = metadata.copy()
        
        if "studies" in metadata:
            # SRA format
            for i, study in enumerate(metadata["studies"]):
                study_context = study.get("study_context", {})
                
                for j, experiment in enumerate(study.get("experiments", [])):
                    seq_type = self.detect_sequencing_type(experiment)
                    extracted_fields = self.extract_sample_fields(
                        experiment, study_context, seq_type
                    )
                    
                    processed_metadata["studies"][i]["experiments"][j]["extracted_metadata"] = extracted_fields
                    processed_metadata["studies"][i]["experiments"][j]["sequencing_type"] = seq_type
                    
        elif "series" in metadata:
            # GEO format
            for i, series in enumerate(metadata["series"]):
                study_context = series.get("study_context", {})
                
                for j, sample in enumerate(series.get("samples", [])):
                    # Infer sequencing type from series metadata
                    seq_type = "RNA-seq"  # Default for GEO
                    extracted_fields = self.extract_sample_fields(
                        sample, study_context, seq_type
                    )
                    
                    processed_metadata["series"][i]["samples"][j]["extracted_metadata"] = extracted_fields
                    processed_metadata["series"][i]["samples"][j]["sequencing_type"] = seq_type
        
        return processed_metadata


def main():
    parser = argparse.ArgumentParser(
        description="Process sample metadata using study context and sequencing-type schemas"
    )
    parser.add_argument("input", help="Input JSON file with contextual metadata")
    parser.add_argument("-o", "--output", help="Output JSON file with processed samples", required=True)
    parser.add_argument("--model-endpoint", help="LLM API endpoint (optional)")
    parser.add_argument("--api-key", help="API key for LLM service (optional)")
    
    args = parser.parse_args()
    
    processor = LLMSampleProcessor(
        model_endpoint=args.model_endpoint,
        api_key=args.api_key
    )
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    processed_metadata = processor.process_samples_with_context(metadata)
    
    with open(args.output, 'w') as f:
        json.dump(processed_metadata, f, indent=2)
    
    print(f"Sample processing completed. Output written to {args.output}")


if __name__ == "__main__":
    main()