#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List

from evidence_system import FieldRegistry, EvidenceCollectionSystem
from evidence_providers import (
    SRAAttributeProvider,
    TitlePatternProvider, 
    GEOCharacteristicsProvider,
    DescriptionPatternProvider
)


class EvidenceBasedMetadataProcessor:
    """Main processor that applies evidence-based standardization to metadata."""
    
    def __init__(self, field_config: str = None):
        # Initialize field registry
        self.field_registry = FieldRegistry(field_config)
        
        # Initialize evidence providers (modular - easy to add new ones)
        self.providers = [
            SRAAttributeProvider(self.field_registry),
            TitlePatternProvider(self.field_registry),
            GEOCharacteristicsProvider(self.field_registry), 
            DescriptionPatternProvider(self.field_registry)
        ]
        
        # Initialize evidence collection system
        self.evidence_system = EvidenceCollectionSystem(
            self.field_registry,
            self.providers
        )
    
    def detect_sequencing_type(self, experiment_metadata: Dict) -> str:
        """Detect sequencing type from experiment metadata (unchanged from original)."""
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
    
    def process_metadata(self, metadata: Dict) -> Dict:
        """Process metadata using evidence-based standardization."""
        processed_metadata = metadata.copy()
        
        # Apply evidence-based standardization to each data unit
        if "studies" in metadata:
            # SRA format - process each experiment
            for i, study in enumerate(metadata["studies"]):
                for j, experiment in enumerate(study.get("experiments", [])):
                    # Detect sequencing type
                    seq_type = self.detect_sequencing_type(experiment)
                    
                    # Apply evidence-based standardization
                    standardized = self.evidence_system.standardize_metadata(experiment)
                    
                    # Add sequencing type and update experiment
                    standardized["sequencing_type"] = seq_type
                    processed_metadata["studies"][i]["experiments"][j] = standardized
                    
        elif "series" in metadata:
            # GEO format - process each sample
            for i, series in enumerate(metadata["series"]):
                for j, sample in enumerate(series.get("samples", [])):
                    # Detect sequencing type (if available)
                    seq_type = self.detect_sequencing_type(sample)
                    
                    # Apply evidence-based standardization 
                    standardized = self.evidence_system.standardize_metadata(sample)
                    
                    # Add sequencing type and update sample
                    standardized["sequencing_type"] = seq_type
                    processed_metadata["series"][i]["samples"][j] = standardized
        
        return processed_metadata


def main():
    parser = argparse.ArgumentParser(
        description="Process metadata using evidence-based standardization"
    )
    parser.add_argument("input", help="Input JSON file with metadata")
    parser.add_argument("-o", "--output", help="Output JSON file with processed metadata", required=True)
    parser.add_argument("--field-config", help="Field registry configuration file")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    
    args = parser.parse_args()
    
    processor = EvidenceBasedMetadataProcessor(args.field_config)
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    if args.debug:
        print("Processing metadata with evidence-based standardization...", file=sys.stderr)
        enabled_fields = processor.field_registry.get_enabled_fields()
        print(f"Enabled fields: {enabled_fields}", file=sys.stderr)
        print(f"Active providers: {[p.get_provider_name() for p in processor.providers]}", file=sys.stderr)
    
    processed_metadata = processor.process_metadata(metadata)
    
    with open(args.output, 'w') as f:
        json.dump(processed_metadata, f, indent=2)
    
    # Print summary statistics
    if "studies" in processed_metadata:
        total_experiments = sum(len(study.get("experiments", [])) for study in processed_metadata["studies"])
        print(f"Processed {total_experiments} experiments using evidence-based standardization")
    elif "series" in processed_metadata:
        total_samples = sum(len(series.get("samples", [])) for series in processed_metadata["series"])
        print(f"Processed {total_samples} samples using evidence-based standardization")


if __name__ == "__main__":
    main()