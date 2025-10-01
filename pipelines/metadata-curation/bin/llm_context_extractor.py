#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Dict, List, Optional, Any
import re


class LLMContextExtractor:
    def __init__(self, model_endpoint: Optional[str] = None, api_key: Optional[str] = None):
        self.model_endpoint = model_endpoint
        self.api_key = api_key
        
    def extract_study_context(self, study_metadata: Dict) -> Dict:
        """Extract study-level experimental context from metadata and publications."""
        context = {
            "study_id": study_metadata.get("bioproject_id") or study_metadata.get("gse_id"),
            "experimental_framework": {},
            "expected_conditions": [],
            "expected_cell_types": [],
            "expected_treatments": [],
            "temporal_design": {},
            "technical_details": {}
        }
        
        # Extract from study description/summary
        description = (study_metadata.get("description") or 
                      study_metadata.get("summary") or 
                      study_metadata.get("overall_design", ""))
        
        if description:
            context["experimental_framework"] = self._parse_experimental_design(description)
            context["expected_conditions"] = self._extract_conditions(description)
            context["expected_cell_types"] = self._extract_cell_types(description)
            context["expected_treatments"] = self._extract_treatments(description)
            context["temporal_design"] = self._extract_temporal_design(description)
            context["technical_details"] = self._extract_technical_details(description)
        
        # Add publication context if available
        if "publication_pmid" in study_metadata or "pubmed_ids" in study_metadata:
            context["publication_context"] = self._extract_publication_context(study_metadata)
        
        return context

    def _parse_experimental_design(self, text: str) -> Dict:
        """Extract experimental design patterns from text."""
        design_patterns = {
            "time_course": bool(re.search(r'\b(time[- ]?course|temporal|longitudinal|over time|time[- ]?point)', text, re.I)),
            "dose_response": bool(re.search(r'\b(dose[- ]?response|concentration|gradient|titration)', text, re.I)),
            "treatment_comparison": bool(re.search(r'\b(treatment|control|versus|compared|treated)', text, re.I)),
            "development": bool(re.search(r'\b(development|developmental|embryonic|stages?)', text, re.I)),
            "differentiation": bool(re.search(r'\b(differentiat|inducd|maturation)', text, re.I)),
            "knockout_knockdown": bool(re.search(r'\b(knockout|knockdown|KO|KD|siRNA|shRNA|CRISPR)', text, re.I)),
            "overexpression": bool(re.search(r'\b(overexpress|ectopic|transfect)', text, re.I))
        }
        
        return design_patterns

    def _extract_conditions(self, text: str) -> List[str]:
        """Extract experimental conditions mentioned in text."""
        conditions = []
        
        # Common condition patterns
        condition_patterns = [
            r'\b(control|untreated|mock|vehicle)\b',
            r'\b(\d+(?:\.\d+)?\s*(?:µM|uM|mM|nM|mg/ml|μg/ml))\b',
            r'\b(treated with|exposed to|stimulated with)\s+([^,\.]+)',
            r'\b(hypoxia|normoxia|stress|starvation)\b',
            r'\b(\d+\s*(?:hours?|hrs?|days?|minutes?|min))\b'
        ]
        
        for pattern in condition_patterns:
            matches = re.findall(pattern, text, re.I)
            conditions.extend([match if isinstance(match, str) else match[1] for match in matches])
        
        return list(set(conditions))

    def _extract_cell_types(self, text: str) -> List[str]:
        """Extract cell types and tissues mentioned in text."""
        cell_types = []
        
        # Common cell type patterns
        cell_patterns = [
            r'\b([A-Z0-9]+\s*cells?)\b',  # HEK293 cells, T cells
            r'\b(fibroblasts?|hepatocytes?|neurons?|cardiomyocytes?)\b',
            r'\b(primary\s+[^,\.]+cells?)\b',
            r'\b([^,\.]+\s+tissue)\b'
        ]
        
        for pattern in cell_patterns:
            matches = re.findall(pattern, text, re.I)
            cell_types.extend(matches)
        
        return list(set(cell_types))

    def _extract_treatments(self, text: str) -> List[str]:
        """Extract treatments and compounds mentioned in text."""
        treatments = []
        
        # Treatment patterns
        treatment_patterns = [
            r'\b([A-Za-z][A-Za-z0-9\-]+)\s*(?:treatment|treated|exposure)',
            r'\b(drug|compound|inhibitor|activator):\s*([^,\.]+)',
            r'\b([A-Za-z]+(?:mycin|cillin|oxin|ide|ase))\b',  # Drug suffixes
        ]
        
        for pattern in treatment_patterns:
            matches = re.findall(pattern, text, re.I)
            treatments.extend([match if isinstance(match, str) else match[1] for match in matches])
        
        return list(set(treatments))

    def _extract_temporal_design(self, text: str) -> Dict:
        """Extract temporal experimental design information."""
        temporal = {
            "is_time_course": False,
            "time_points": [],
            "duration": None
        }
        
        # Time point patterns
        time_patterns = [
            r'\b(\d+(?:\.\d+)?)\s*(hours?|hrs?|h)\b',
            r'\b(\d+(?:\.\d+)?)\s*(days?|d)\b',
            r'\b(\d+(?:\.\d+)?)\s*(weeks?|w)\b',
            r'\b(\d+(?:\.\d+)?)\s*(minutes?|min|m)\b'
        ]
        
        for pattern in time_patterns:
            matches = re.findall(pattern, text, re.I)
            for value, unit in matches:
                temporal["time_points"].append(f"{value} {unit}")
                temporal["is_time_course"] = True
        
        return temporal

    def _extract_technical_details(self, text: str) -> Dict:
        """Extract technical sequencing details."""
        technical = {
            "library_type": None,
            "sequencing_platform": None,
            "paired_end": None,
            "read_length": None,
            "depth": None
        }
        
        # Library type
        lib_patterns = [
            r'\b(RNA[- ]?seq|ribosome profiling|ribo[- ]?seq|ChIP[- ]?seq|ATAC[- ]?seq)\b'
        ]
        for pattern in lib_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                technical["library_type"] = match.group(1)
                break
        
        # Platform
        platform_patterns = [
            r'\b(Illumina|HiSeq|NextSeq|NovaSeq|MiSeq|PacBio|Oxford Nanopore|Ion Torrent)\b'
        ]
        for pattern in platform_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                technical["sequencing_platform"] = match.group(1)
                break
        
        # Read configuration
        if re.search(r'\b(paired[- ]?end|PE)\b', text, re.I):
            technical["paired_end"] = True
        elif re.search(r'\b(single[- ]?end|SE)\b', text, re.I):
            technical["paired_end"] = False
            
        # Read length
        read_match = re.search(r'\b(\d+)\s*(?:bp|nt|nucleotides?|base[- ]?pairs?)\b', text, re.I)
        if read_match:
            technical["read_length"] = int(read_match.group(1))
        
        return technical

    def _extract_publication_context(self, study_metadata: Dict) -> Dict:
        """Extract context from linked publications."""
        pub_context = {
            "pmid": study_metadata.get("publication_pmid"),
            "has_methods": False,
            "experimental_description": None
        }
        
        # This would integrate with PMC extraction results
        # For now, return basic structure
        return pub_context

    def process_metadata_with_context(self, metadata: Dict) -> Dict:
        """Process metadata and add contextual information."""
        enhanced_metadata = metadata.copy()
        
        if "studies" in metadata:
            # SRA format
            for i, study in enumerate(metadata["studies"]):
                bioproject_meta = study.get("bioproject_metadata", {})
                context = self.extract_study_context(bioproject_meta)
                enhanced_metadata["studies"][i]["study_context"] = context
                
        elif "series" in metadata:
            # GEO format
            for i, series in enumerate(metadata["series"]):
                series_meta = series.get("series_metadata", {})
                context = self.extract_study_context(series_meta)
                enhanced_metadata["series"][i]["study_context"] = context
                
        elif "articles" in metadata:
            # PMC format - extract context from full text
            for i, article in enumerate(metadata["articles"]):
                context = self._extract_article_context(article)
                enhanced_metadata["articles"][i]["extracted_context"] = context
        
        return enhanced_metadata

    def _extract_article_context(self, article: Dict) -> Dict:
        """Extract experimental context from article content."""
        context = {
            "experimental_methods": [],
            "sample_processing": [],
            "sequencing_protocols": [],
            "data_analysis": []
        }
        
        methods_text = article.get("methods", "")
        if methods_text:
            context["experimental_methods"] = self._parse_methods_text(methods_text)
            
        materials_text = article.get("materials", "")
        if materials_text:
            context["sample_processing"] = self._parse_sample_processing(materials_text)
            
        return context

    def _parse_methods_text(self, methods_text: str) -> List[str]:
        """Parse methods section for key experimental details."""
        methods = []
        
        # Extract sentence-level methods
        sentences = re.split(r'[.!?]+', methods_text)
        for sentence in sentences:
            sentence = sentence.strip()
            if len(sentence) > 20:  # Filter out very short fragments
                # Look for method-related keywords
                if re.search(r'\b(cultured|incubated|treated|extracted|sequenced|prepared)\b', sentence, re.I):
                    methods.append(sentence)
        
        return methods[:10]  # Limit to top 10 methods

    def _parse_sample_processing(self, materials_text: str) -> List[str]:
        """Parse sample processing details."""
        processing = []
        
        sentences = re.split(r'[.!?]+', materials_text)
        for sentence in sentences:
            sentence = sentence.strip()
            if len(sentence) > 20:
                if re.search(r'\b(RNA|library|protocol|extraction|purification)\b', sentence, re.I):
                    processing.append(sentence)
        
        return processing[:10]  # Limit to top 10 processing steps


def main():
    parser = argparse.ArgumentParser(
        description="Extract experimental context from study metadata using pattern matching"
    )
    parser.add_argument("input", help="Input JSON file with metadata")
    parser.add_argument("-o", "--output", help="Output JSON file with context", required=True)
    parser.add_argument("--model-endpoint", help="LLM API endpoint (optional)")
    parser.add_argument("--api-key", help="API key for LLM service (optional)")
    
    args = parser.parse_args()
    
    extractor = LLMContextExtractor(
        model_endpoint=args.model_endpoint,
        api_key=args.api_key
    )
    
    with open(args.input) as f:
        metadata = json.load(f)
    
    enhanced_metadata = extractor.process_metadata_with_context(metadata)
    
    with open(args.output, 'w') as f:
        json.dump(enhanced_metadata, f, indent=2)
    
    print(f"Context extraction completed. Output written to {args.output}")


if __name__ == "__main__":
    main()