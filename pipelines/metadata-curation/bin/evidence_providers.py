#!/usr/bin/env python3

import re
from typing import Dict, List, Set
from evidence_system import EvidenceProvider, Evidence, FieldRegistry


class SRAAttributeProvider(EvidenceProvider):
    """Extract evidence from SRA sample_attributes fields."""
    
    def __init__(self, field_registry: FieldRegistry):
        self.field_registry = field_registry
    
    def get_provider_name(self) -> str:
        return "SRAAttributeProvider"
    
    def get_supported_fields(self) -> Set[str]:
        return {"tissue", "cell_type", "treatment", "timepoint", "genotype", "replicate"}
    
    def extract_evidence(self, metadata: Dict, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        # Check if this is a single experiment (from processing pipeline)
        if "sample_attributes" in metadata:
            sample_attrs = metadata["sample_attributes"]
            evidence.extend(self._extract_from_attributes(sample_attrs, target_fields))
        
        # Also handle full study format (from raw data)
        elif "studies" in metadata:
            for study in metadata["studies"]:
                for experiment in study.get("experiments", []):
                    sample_attrs = experiment.get("sample_attributes", {})
                    evidence.extend(self._extract_from_attributes(sample_attrs, target_fields))
        
        return evidence
    
    def _extract_from_attributes(self, attributes: Dict, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        # Direct field name matching
        for attr_name, attr_value in attributes.items():
            if not isinstance(attr_value, str) or not attr_value.strip():
                continue
                
            attr_lower = attr_name.lower()
            
            # Check each target field's search patterns
            for target_field in target_fields:
                patterns = self.field_registry.get_search_patterns(target_field)
                
                for pattern in patterns:
                    if pattern in attr_lower:
                        evidence.append(Evidence(
                            field_target=target_field,
                            source_type="sample_attribute",
                            source_field=attr_name,
                            raw_value=attr_value,
                            extraction_method="direct_field_match",
                            provider=self.get_provider_name()
                        ))
                        break  # Only one match per field per attribute
        
        return evidence


class TitlePatternProvider(EvidenceProvider):
    """Extract evidence from experiment and study titles using patterns."""
    
    def __init__(self, field_registry: FieldRegistry):
        self.field_registry = field_registry
    
    def get_provider_name(self) -> str:
        return "TitlePatternProvider"
    
    def get_supported_fields(self) -> Set[str]:
        return {"tissue", "cell_type", "treatment", "timepoint"}
    
    def extract_evidence(self, metadata: Dict, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        # Extract from SRA experiment titles
        if "studies" in metadata:
            for study in metadata["studies"]:
                study_title = study.get("bioproject_metadata", {}).get("title", "")
                if study_title:
                    evidence.extend(self._extract_from_title(study_title, "study_title", target_fields))
                
                for experiment in study.get("experiments", []):
                    exp_title = experiment.get("experiment_title", "")
                    if exp_title:
                        evidence.extend(self._extract_from_title(exp_title, "experiment_title", target_fields))
        
        # Extract from GEO series titles
        if "series" in metadata:
            for series in metadata["series"]:
                series_title = series.get("series_metadata", {}).get("title", "")
                if series_title:
                    evidence.extend(self._extract_from_title(series_title, "series_title", target_fields))
        
        return evidence
    
    def _extract_from_title(self, title: str, source_field: str, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        # Tissue extraction patterns
        if "tissue" in target_fields:
            tissue_patterns = [
                r'(\w+)\s+tissue',
                r'(\w+)\s+cells?',
                r'in\s+(\w+)',
                r'from\s+(\w+)',
                r'(\w+)\s+samples?'
            ]
            
            for pattern in tissue_patterns:
                matches = re.findall(pattern, title.lower())
                for match in matches:
                    if len(match) > 2:  # Skip very short matches
                        evidence.append(Evidence(
                            field_target="tissue",
                            source_type="title",
                            source_field=source_field,
                            raw_value=match,
                            extraction_method="regex_pattern",
                            provider=self.get_provider_name()
                        ))
        
        # Cell type extraction patterns  
        if "cell_type" in target_fields:
            cell_patterns = [
                r'(hek\d+[a-z]*)',
                r'(hela)',
                r'(cho)',
                r'(\w+)\s+cell\s+line',
                r'(\w+)\s+cells?'
            ]
            
            for pattern in cell_patterns:
                matches = re.findall(pattern, title.lower())
                for match in matches:
                    evidence.append(Evidence(
                        field_target="cell_type",
                        source_type="title",
                        source_field=source_field,
                        raw_value=match,
                        extraction_method="regex_pattern",
                        provider=self.get_provider_name()
                    ))
        
        # Treatment extraction patterns
        if "treatment" in target_fields:
            treatment_patterns = [
                r'treated\s+with\s+(\w+)',
                r'(\w+)\s+treatment',
                r'plus\s+(\w+)',
                r'(\w+)\s+treated'
            ]
            
            for pattern in treatment_patterns:
                matches = re.findall(pattern, title.lower())
                for match in matches:
                    evidence.append(Evidence(
                        field_target="treatment",
                        source_type="title", 
                        source_field=source_field,
                        raw_value=match,
                        extraction_method="regex_pattern",
                        provider=self.get_provider_name()
                    ))
        
        return evidence


class GEOCharacteristicsProvider(EvidenceProvider):
    """Extract evidence from GEO characteristics fields."""
    
    def __init__(self, field_registry: FieldRegistry):
        self.field_registry = field_registry
    
    def get_provider_name(self) -> str:
        return "GEOCharacteristicsProvider"
    
    def get_supported_fields(self) -> Set[str]:
        return {"tissue", "cell_type", "treatment", "timepoint", "genotype"}
    
    def extract_evidence(self, metadata: Dict, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        if "series" in metadata:
            for series in metadata["series"]:
                for sample in series.get("samples", []):
                    characteristics = sample.get("characteristics", {})
                    evidence.extend(self._extract_from_characteristics(characteristics, target_fields))
        
        return evidence
    
    def _extract_from_characteristics(self, characteristics: Dict, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        for char_name, char_value in characteristics.items():
            if not isinstance(char_value, str) or not char_value.strip():
                continue
            
            char_lower = char_name.lower()
            
            # Check each target field's search patterns
            for target_field in target_fields:
                patterns = self.field_registry.get_search_patterns(target_field)
                
                for pattern in patterns:
                    if pattern in char_lower:
                        evidence.append(Evidence(
                            field_target=target_field,
                            source_type="geo_characteristic",
                            source_field=char_name,
                            raw_value=char_value,
                            extraction_method="direct_field_match",
                            provider=self.get_provider_name()
                        ))
                        break
        
        return evidence


class DescriptionPatternProvider(EvidenceProvider):
    """Extract evidence from free text descriptions using NLP patterns."""
    
    def __init__(self, field_registry: FieldRegistry):
        self.field_registry = field_registry
    
    def get_provider_name(self) -> str:
        return "DescriptionPatternProvider"
    
    def get_supported_fields(self) -> Set[str]:
        return {"tissue", "cell_type", "treatment", "timepoint"}
    
    def extract_evidence(self, metadata: Dict, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        # Extract from study descriptions
        if "studies" in metadata:
            for study in metadata["studies"]:
                description = study.get("bioproject_metadata", {}).get("description", "")
                if description:
                    evidence.extend(self._extract_from_description(description, "study_description", target_fields))
        
        return evidence
    
    def _extract_from_description(self, description: str, source_field: str, target_fields: Set[str]) -> List[Evidence]:
        evidence = []
        
        # Simple pattern-based extraction (could be enhanced with real NLP)
        desc_lower = description.lower()
        
        # Tissue patterns in descriptions
        if "tissue" in target_fields:
            tissue_patterns = [
                r'(\w+)\s+tissue\s+samples?',
                r'samples?\s+from\s+(\w+)',
                r'(\w+)\s+biopsy',
                r'(\w+)\s+organ'
            ]
            
            for pattern in tissue_patterns:
                matches = re.findall(pattern, desc_lower)
                for match in matches:
                    if len(match) > 2:
                        evidence.append(Evidence(
                            field_target="tissue",
                            source_type="description",
                            source_field=source_field,
                            raw_value=match,
                            extraction_method="description_pattern",
                            provider=self.get_provider_name()
                        ))
        
        return evidence