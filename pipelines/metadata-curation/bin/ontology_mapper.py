#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from pathlib import Path
import re

import requests


@dataclass
class OntologyTerm:
    id: str
    name: str
    synonyms: List[str]
    definition: Optional[str] = None
    
    
@dataclass
class MappingResult:
    original_term: str
    mapped_term: Optional[OntologyTerm]
    confidence: str  # "exact", "synonym", "unmapped"
    source_ontology: Optional[str] = None


class OntologyMapper:
    def __init__(self, cache_dir: str = "ontology_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Core ontologies for biomedical metadata
        self.ontologies = {
            "CHEBI": "https://www.ebi.ac.uk/ols4/api/ontologies/chebi",
            "CL": "https://www.ebi.ac.uk/ols4/api/ontologies/cl", 
            "EFO": "https://www.ebi.ac.uk/ols4/api/ontologies/efo",
            "NCBITAXON": "https://www.ebi.ac.uk/ols4/api/ontologies/ncbitaxon",
            "GO": "https://www.ebi.ac.uk/ols4/api/ontologies/go",
            "HP": "https://www.ebi.ac.uk/ols4/api/ontologies/hp",
            "UBERON": "https://www.ebi.ac.uk/ols4/api/ontologies/uberon"
        }
        
        # Preloaded synonym mappings for common terms
        self.manual_synonyms = {
            "cycloheximide": ["cyclohexamide", "CHX", "actidione"],
            "dimethyl sulfoxide": ["DMSO", "dimethylsulfoxide"],
            "fetal bovine serum": ["FBS", "foetal bovine serum", "fetal calf serum"],
            "phosphate buffered saline": ["PBS", "phosphate-buffered saline"],
            "trizol": ["TRIzol", "tri-reagent"],
            "untreated": ["control", "mock", "vehicle", "no treatment"]
        }

    def load_ontology_terms(self, ontology: str, max_terms: int = 10000) -> List[OntologyTerm]:
        """Load terms from an ontology via OLS API."""
        cache_file = self.cache_dir / f"{ontology.lower()}_terms.json"
        
        if cache_file.exists():
            with open(cache_file) as f:
                cached_data = json.load(f)
                return [OntologyTerm(**term) for term in cached_data]
        
        print(f"Loading {ontology} ontology terms...", file=sys.stderr)
        terms = []
        
        try:
            url = f"{self.ontologies[ontology]}/terms"
            params = {"size": max_terms}
            
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            data = response.json()
            
            for term_data in data.get("_embedded", {}).get("terms", []):
                term = OntologyTerm(
                    id=term_data.get("obo_id", ""),
                    name=term_data.get("label", ""),
                    synonyms=term_data.get("synonyms", []),
                    definition=term_data.get("description", [""])[0] if term_data.get("description") else None
                )
                terms.append(term)
            
            # Cache the results
            with open(cache_file, 'w') as f:
                json.dump([term.__dict__ for term in terms], f)
                
        except Exception as e:
            print(f"Error loading {ontology}: {e}", file=sys.stderr)
            
        return terms

    def normalize_term(self, term: str) -> str:
        """Normalize a term for matching."""
        if not term:
            return ""
        
        # Convert to lowercase and clean
        normalized = term.lower().strip()
        
        # Remove common prefixes/suffixes
        normalized = re.sub(r'^(treatment|condition):\s*', '', normalized)
        normalized = re.sub(r'\s*(treatment|condition)$', '', normalized)
        
        # Standardize units and formatting
        normalized = re.sub(r'\s+', ' ', normalized)  # Multiple spaces to single
        normalized = re.sub(r'[^\w\s\-\.]', '', normalized)  # Remove special chars except dash, dot
        
        return normalized

    def exact_match(self, term: str, ontology_terms: List[OntologyTerm]) -> Optional[OntologyTerm]:
        """Find exact match in ontology terms or synonyms."""
        normalized_term = self.normalize_term(term)
        
        for ont_term in ontology_terms:
            # Check primary label
            if self.normalize_term(ont_term.name) == normalized_term:
                return ont_term
                
            # Check synonyms
            for synonym in ont_term.synonyms:
                if self.normalize_term(synonym) == normalized_term:
                    return ont_term
                    
        return None

    def manual_synonym_match(self, term: str) -> Optional[str]:
        """Check against manually curated synonym mappings."""
        normalized_term = self.normalize_term(term)
        
        for canonical_term, synonyms in self.manual_synonyms.items():
            if normalized_term == self.normalize_term(canonical_term):
                return canonical_term
            for synonym in synonyms:
                if normalized_term == self.normalize_term(synonym):
                    return canonical_term
                    
        return None

    def map_term(self, term: str, ontology_terms: Dict[str, List[OntologyTerm]]) -> MappingResult:
        """Map a single term to ontology concepts."""
        if not term or not term.strip():
            return MappingResult(term, None, "unmapped")
        
        # First check manual synonyms
        manual_match = self.manual_synonym_match(term)
        if manual_match:
            # Try to find the canonical term in ontologies
            for ont_name, ont_terms in ontology_terms.items():
                exact_term = self.exact_match(manual_match, ont_terms)
                if exact_term:
                    return MappingResult(term, exact_term, "synonym", ont_name)
        
        # Try exact matching across all ontologies
        for ont_name, ont_terms in ontology_terms.items():
            exact_term = self.exact_match(term, ont_terms)
            if exact_term:
                return MappingResult(term, exact_term, "exact", ont_name)
        
        return MappingResult(term, None, "unmapped")

    def map_metadata_fields(self, metadata: Dict, field_mappings: Dict[str, str]) -> Dict:
        """Map metadata fields to ontology terms."""
        # Load all required ontologies
        ontology_terms = {}
        for ont_name in self.ontologies.keys():
            ontology_terms[ont_name] = self.load_ontology_terms(ont_name)
        
        mapped_metadata = metadata.copy()
        mapping_results = {}
        
        for field_name, ont_target in field_mappings.items():
            if field_name in metadata:
                field_value = metadata[field_name]
                
                if isinstance(field_value, str):
                    result = self.map_term(field_value, ontology_terms)
                    mapping_results[field_name] = result
                    
                    if result.mapped_term:
                        mapped_metadata[f"{field_name}_ontology_id"] = result.mapped_term.id
                        mapped_metadata[f"{field_name}_ontology_name"] = result.mapped_term.name
                        mapped_metadata[f"{field_name}_ontology_source"] = result.source_ontology
                    
                elif isinstance(field_value, dict):
                    # Handle nested dictionaries (like sample_attributes)
                    for key, value in field_value.items():
                        if isinstance(value, str):
                            result = self.map_term(value, ontology_terms)
                            mapping_results[f"{field_name}.{key}"] = result
                            
                            if result.mapped_term:
                                if f"{field_name}_mapped" not in mapped_metadata:
                                    mapped_metadata[f"{field_name}_mapped"] = {}
                                mapped_metadata[f"{field_name}_mapped"][key] = {
                                    "ontology_id": result.mapped_term.id,
                                    "ontology_name": result.mapped_term.name,
                                    "ontology_source": result.source_ontology,
                                    "confidence": result.confidence
                                }
        
        mapped_metadata["_mapping_results"] = {
            field: {
                "original": result.original_term,
                "mapped_id": result.mapped_term.id if result.mapped_term else None,
                "mapped_name": result.mapped_term.name if result.mapped_term else None,
                "confidence": result.confidence,
                "source_ontology": result.source_ontology
            }
            for field, result in mapping_results.items()
        }
        
        return mapped_metadata

    def generate_mapping_report(self, mapping_results: Dict[str, MappingResult]) -> Dict:
        """Generate a summary report of mapping quality."""
        total_terms = len(mapping_results)
        exact_matches = sum(1 for r in mapping_results.values() if r.confidence == "exact")
        synonym_matches = sum(1 for r in mapping_results.values() if r.confidence == "synonym")
        unmapped = sum(1 for r in mapping_results.values() if r.confidence == "unmapped")
        
        unmapped_terms = [
            result.original_term 
            for result in mapping_results.values() 
            if result.confidence == "unmapped"
        ]
        
        # Count frequency of unmapped terms
        unmapped_freq = {}
        for term in unmapped_terms:
            unmapped_freq[term] = unmapped_freq.get(term, 0) + 1
        
        sorted_unmapped = sorted(unmapped_freq.items(), key=lambda x: x[1], reverse=True)
        
        return {
            "total_terms": total_terms,
            "exact_matches": exact_matches,
            "synonym_matches": synonym_matches,
            "unmapped": unmapped,
            "mapping_rate": (exact_matches + synonym_matches) / total_terms if total_terms > 0 else 0,
            "top_unmapped_terms": sorted_unmapped[:20]  # Top 20 most frequent unmapped terms
        }


def main():
    parser = argparse.ArgumentParser(
        description="Map metadata terms to ontology concepts"
    )
    parser.add_argument("input", help="Input JSON file with metadata")
    parser.add_argument("-o", "--output", help="Output JSON file with mapped metadata", required=True)
    parser.add_argument("--field-mappings", help="JSON file specifying field->ontology mappings")
    parser.add_argument("--cache-dir", default="ontology_cache", 
                        help="Directory for ontology cache files")
    
    args = parser.parse_args()
    
    # Default field mappings for common metadata fields
    default_mappings = {
        "organism": "NCBITAXON",
        "cell_type": "CL",
        "tissue": "UBERON",
        "disease": "HP",
        "treatment": "CHEBI",
        "cell_line": "EFO"
    }
    
    if args.field_mappings:
        with open(args.field_mappings) as f:
            field_mappings = json.load(f)
    else:
        field_mappings = default_mappings
    
    mapper = OntologyMapper(cache_dir=args.cache_dir)
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    # Map metadata recursively
    if isinstance(metadata, dict):
        if "studies" in metadata:
            # SRA format
            for study in metadata["studies"]:
                for experiment in study.get("experiments", []):
                    experiment.update(mapper.map_metadata_fields(experiment, field_mappings))
        elif "series" in metadata:
            # GEO format  
            for series in metadata["series"]:
                for sample in series.get("samples", []):
                    sample.update(mapper.map_metadata_fields(sample, field_mappings))
        else:
            # Single metadata object
            metadata.update(mapper.map_metadata_fields(metadata, field_mappings))
    
    with open(args.output, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"Ontology mapping completed. Output written to {args.output}")


if __name__ == "__main__":
    main()