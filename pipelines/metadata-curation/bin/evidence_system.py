#!/usr/bin/env python3

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Any
import re
import json


@dataclass
class Evidence:
    """Single piece of evidence for a standardized field."""
    field_target: str       # What standardized field this relates to (e.g., "tissue")
    source_type: str        # Where it came from (e.g., "sample_attribute") 
    source_field: str       # Exact field name (e.g., "tissue_type")
    raw_value: str          # Original value (e.g., "kidney cortex")
    extraction_method: str  # How it was found (e.g., "direct_field_match")
    provider: str           # Which provider found it


@dataclass
class StandardizedField:
    """Result of evidence resolution for a standardized field."""
    field_name: str         # Standardized field name
    standardized_value: str # Final standardized value
    evidence_used: Evidence # Primary evidence that determined the value
    all_evidence: List[Evidence]  # All evidence collected
    rejected_evidence: List[Evidence]  # Evidence that was rejected
    resolution_method: str  # How conflicts were resolved


class EvidenceProvider(ABC):
    """Abstract base class for evidence providers."""
    
    @abstractmethod
    def get_provider_name(self) -> str:
        """Return name of this provider."""
        pass
    
    @abstractmethod
    def get_supported_fields(self) -> Set[str]:
        """Return set of standardized fields this provider can find evidence for."""
        pass
    
    @abstractmethod
    def extract_evidence(self, metadata: Dict, target_fields: Set[str]) -> List[Evidence]:
        """Extract evidence for target fields from metadata."""
        pass


class FieldRegistry:
    """Registry of standardized fields and their search patterns."""
    
    def __init__(self, config_file: Optional[str] = None):
        # Default field definitions
        self.field_definitions = {
            "tissue": {
                "description": "Anatomical tissue or organ",
                "ontology": "UBERON",
                "search_patterns": ["tissue", "organ", "source_tissue", "tissue_type"],
                "enabled": True
            },
            "cell_type": {
                "description": "Cell type or cell line", 
                "ontology": "CL",
                "search_patterns": ["cell_type", "cell_line", "cell", "cells"],
                "enabled": True
            },
            "treatment": {
                "description": "Chemical treatment or compound",
                "ontology": "CHEBI", 
                "search_patterns": ["treatment", "drug", "compound", "chemical", "reagent"],
                "enabled": True
            },
            "timepoint": {
                "description": "Time point of sample collection",
                "ontology": None,
                "search_patterns": ["time", "timepoint", "time_point", "duration", "hours", "days"],
                "enabled": True
            },
            "genotype": {
                "description": "Genetic background or strain",
                "ontology": None,
                "search_patterns": ["genotype", "strain", "background", "genetic_background"],
                "enabled": True
            },
            "replicate": {
                "description": "Biological or technical replicate",
                "ontology": None,
                "search_patterns": ["replicate", "rep", "biological_replicate", "technical_replicate"],
                "enabled": True
            }
        }
        
        if config_file:
            self.load_config(config_file)
    
    def load_config(self, config_file: str):
        """Load field definitions from config file."""
        with open(config_file) as f:
            config = json.load(f)
            self.field_definitions.update(config.get("fields", {}))
    
    def get_enabled_fields(self) -> Set[str]:
        """Get set of enabled standardized field names."""
        return {name for name, defn in self.field_definitions.items() 
                if defn.get("enabled", True)}
    
    def get_search_patterns(self, field_name: str) -> List[str]:
        """Get search patterns for a standardized field."""
        return self.field_definitions.get(field_name, {}).get("search_patterns", [])
    
    def get_ontology(self, field_name: str) -> Optional[str]:
        """Get target ontology for a standardized field."""
        return self.field_definitions.get(field_name, {}).get("ontology")


class EvidenceResolver:
    """Resolves conflicts between evidence to produce standardized fields."""
    
    def __init__(self, field_registry: FieldRegistry):
        self.field_registry = field_registry
        
        # Define resolution rules (transparent, no black box scoring)
        self.source_priority = [
            "sample_attribute",
            "experiment_attribute", 
            "experiment_title",
            "study_title",
            "study_description"
        ]
        
        self.extraction_priority = [
            "direct_field_match",
            "pattern_match",
            "keyword_extraction",
            "fuzzy_match"
        ]
    
    def resolve_field_evidence(self, field_name: str, evidence_list: List[Evidence]) -> StandardizedField:
        """Resolve evidence for a single standardized field."""
        if not evidence_list:
            return StandardizedField(
                field_name=field_name,
                standardized_value="unknown",
                evidence_used=None,
                all_evidence=[],
                rejected_evidence=[],
                resolution_method="no_evidence"
            )
        
        # Sort evidence by priority rules
        sorted_evidence = self._sort_evidence_by_priority(evidence_list)
        
        # Take the highest priority evidence
        primary_evidence = sorted_evidence[0]
        rejected = sorted_evidence[1:] if len(sorted_evidence) > 1 else []
        
        # Apply any value standardization (future: use ontology mapping here)
        standardized_value = self._standardize_value(field_name, primary_evidence.raw_value)
        
        return StandardizedField(
            field_name=field_name,
            standardized_value=standardized_value,
            evidence_used=primary_evidence,
            all_evidence=evidence_list,
            rejected_evidence=rejected,
            resolution_method=f"priority_rule_{primary_evidence.source_type}_{primary_evidence.extraction_method}"
        )
    
    def _sort_evidence_by_priority(self, evidence_list: List[Evidence]) -> List[Evidence]:
        """Sort evidence by priority rules."""
        def priority_score(evidence: Evidence) -> tuple:
            source_idx = self.source_priority.index(evidence.source_type) if evidence.source_type in self.source_priority else 999
            extraction_idx = self.extraction_priority.index(evidence.extraction_method) if evidence.extraction_method in self.extraction_priority else 999
            return (source_idx, extraction_idx)
        
        return sorted(evidence_list, key=priority_score)
    
    def _standardize_value(self, field_name: str, raw_value: str) -> str:
        """Apply basic value standardization (placeholder for ontology integration)."""
        if not raw_value or raw_value.lower() in ["unknown", "not applicable", "n/a", ""]:
            return "unknown"
        
        # Basic cleanup
        cleaned = raw_value.strip().lower()
        
        # Field-specific standardization rules (expandable)
        if field_name == "timepoint":
            # Standardize time formats
            cleaned = re.sub(r'(\d+)\s*hrs?', r'\1 hours', cleaned)
            cleaned = re.sub(r'(\d+)\s*h\b', r'\1 hours', cleaned)
            cleaned = re.sub(r'(\d+)\s*mins?', r'\1 minutes', cleaned)
            cleaned = re.sub(r'(\d+)\s*m\b', r'\1 minutes', cleaned)
        
        return cleaned


class EvidenceCollectionSystem:
    """Main system that orchestrates evidence collection and resolution."""
    
    def __init__(self, field_registry: FieldRegistry, providers: List[EvidenceProvider]):
        self.field_registry = field_registry
        self.providers = providers
        self.resolver = EvidenceResolver(field_registry)
    
    def standardize_metadata(self, metadata: Dict) -> Dict:
        """Apply evidence-based standardization to metadata."""
        target_fields = self.field_registry.get_enabled_fields()
        
        # Collect evidence from all providers
        all_evidence = []
        for provider in self.providers:
            evidence = provider.extract_evidence(metadata, target_fields)
            all_evidence.extend(evidence)
        
        # Group evidence by target field
        evidence_by_field = {}
        for evidence in all_evidence:
            if evidence.field_target not in evidence_by_field:
                evidence_by_field[evidence.field_target] = []
            evidence_by_field[evidence.field_target].append(evidence)
        
        # Resolve evidence for each field
        standardized_fields = {}
        for field_name in target_fields:
            field_evidence = evidence_by_field.get(field_name, [])
            resolved = self.resolver.resolve_field_evidence(field_name, field_evidence)
            standardized_fields[field_name] = resolved
        
        # Create output with standardized metadata
        result = metadata.copy()
        result["_evidence_standardization"] = {
            "standardized_fields": {
                field_name: {
                    "value": resolved.standardized_value,
                    "evidence_used": {
                        "source": f"{resolved.evidence_used.source_type}.{resolved.evidence_used.source_field}",
                        "raw_value": resolved.evidence_used.raw_value,
                        "method": resolved.evidence_used.extraction_method,
                        "provider": resolved.evidence_used.provider
                    } if resolved.evidence_used else None,
                    "all_evidence_count": len(resolved.all_evidence),
                    "rejected_evidence_count": len(resolved.rejected_evidence),
                    "resolution_method": resolved.resolution_method
                }
                for field_name, resolved in standardized_fields.items()
            }
        }
        
        # Add standardized values to main metadata
        for field_name, resolved in standardized_fields.items():
            result[f"standardized_{field_name}"] = resolved.standardized_value
        
        return result