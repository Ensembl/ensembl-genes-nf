#!/usr/bin/env python3

import argparse
import json
import sys
import time
from typing import Dict, List, Optional, Set
from xml.etree import ElementTree as ET

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class SRAMetadataExtractor:
    def __init__(self, email: Optional[str] = None, max_retries: int = 3):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.email = email
        self.session = self._setup_session(max_retries)
        
    def _setup_session(self, max_retries: int) -> requests.Session:
        session = requests.Session()
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        return session

    def search_studies(self, query: str, max_results: int = 1000) -> List[str]:
        """Search for SRA studies matching the query and return BioProject IDs."""
        params = {
            "db": "bioproject",
            "term": f"{query}[Description] OR {query}[Title]",
            "retmax": max_results,
            "retmode": "xml"
        }
        if self.email:
            params["email"] = self.email
            
        response = self._make_request("esearch.fcgi", params)
        root = ET.fromstring(response.content)
        
        ids = []
        for id_elem in root.findall(".//Id"):
            ids.append(id_elem.text)
            
        return ids

    def get_study_metadata(self, bioproject_id: str) -> Dict:
        """Fetch detailed metadata for a BioProject."""
        params = {
            "db": "bioproject",
            "id": bioproject_id,
            "retmode": "xml"
        }
        if self.email:
            params["email"] = self.email
            
        response = self._make_request("efetch.fcgi", params)
        root = ET.fromstring(response.content)
        
        project_elem = root.find(".//Project")
        if project_elem is None:
            return {}
            
        metadata = {
            "bioproject_id": bioproject_id,
            "title": self._get_text(project_elem, ".//Title"),
            "description": self._get_text(project_elem, ".//Description"),
            "organism": self._get_text(project_elem, ".//OrganismName"),
            "taxonomy_id": self._get_text(project_elem, ".//TaxonomyId"),
            "data_type": self._get_text(project_elem, ".//ProjectDataType"),
            "publication_pmid": self._get_text(project_elem, ".//PublicationId[@id_type='pubmed']"),
            "submission_date": self._get_text(project_elem, ".//SubmissionDate"),
            "release_date": self._get_text(project_elem, ".//ReleaseDate")
        }
        
        return metadata

    def get_linked_experiments(self, bioproject_id: str) -> List[str]:
        """Get SRA experiment IDs linked to a BioProject."""
        params = {
            "dbfrom": "bioproject",
            "db": "sra",
            "id": bioproject_id,
            "retmode": "xml"
        }
        if self.email:
            params["email"] = self.email
            
        response = self._make_request("elink.fcgi", params)
        root = ET.fromstring(response.content)
        
        experiment_ids = []
        for id_elem in root.findall(".//Id"):
            experiment_ids.append(id_elem.text)
            
        return experiment_ids

    def get_experiment_metadata(self, experiment_ids: List[str]) -> List[Dict]:
        """Fetch detailed metadata for SRA experiments."""
        if not experiment_ids:
            return []
            
        params = {
            "db": "sra",
            "id": ",".join(experiment_ids),
            "retmode": "xml"
        }
        if self.email:
            params["email"] = self.email
            
        response = self._make_request("efetch.fcgi", params)
        root = ET.fromstring(response.content)
        
        experiments = []
        for exp_elem in root.findall(".//EXPERIMENT_PACKAGE"):
            metadata = self._parse_experiment(exp_elem)
            if metadata:
                experiments.append(metadata)
                
        return experiments

    def _parse_experiment(self, exp_elem: ET.Element) -> Optional[Dict]:
        """Parse experiment XML element into metadata dictionary."""
        try:
            exp = exp_elem.find(".//EXPERIMENT")
            sample = exp_elem.find(".//SAMPLE")
            study = exp_elem.find(".//STUDY")
            run_set = exp_elem.find(".//RUN_SET")
            
            if exp is None:
                return None
                
            metadata = {
                "experiment_id": exp.get("accession"),
                "experiment_title": self._get_text(exp, ".//TITLE"),
                "study_id": study.get("accession") if study is not None else None,
                "study_title": self._get_text(study, ".//STUDY_TITLE"),
                "sample_id": sample.get("accession") if sample is not None else None,
                "sample_title": self._get_text(sample, ".//TITLE"),
                "organism": self._get_text(sample, ".//SCIENTIFIC_NAME"),
                "library_strategy": self._get_text(exp, ".//LIBRARY_STRATEGY"),
                "library_source": self._get_text(exp, ".//LIBRARY_SOURCE"),
                "library_selection": self._get_text(exp, ".//LIBRARY_SELECTION"),
                "platform": self._get_text(exp, ".//PLATFORM_NAME"),
                "instrument": self._get_text(exp, ".//INSTRUMENT_MODEL"),
                "library_layout": "paired" if exp.find(".//PAIRED") is not None else "single",
                "runs": []
            }
            
            # Parse sample attributes
            sample_attrs = {}
            if sample is not None:
                for attr in sample.findall(".//SAMPLE_ATTRIBUTE"):
                    tag = self._get_text(attr, ".//TAG")
                    value = self._get_text(attr, ".//VALUE")
                    if tag and value:
                        sample_attrs[tag] = value
            metadata["sample_attributes"] = sample_attrs
            
            # Parse runs
            if run_set is not None:
                for run in run_set.findall(".//RUN"):
                    run_metadata = {
                        "run_id": run.get("accession"),
                        "total_spots": run.get("total_spots"),
                        "total_bases": run.get("total_bases"),
                        "size": run.get("size"),
                        "published": run.get("published")
                    }
                    metadata["runs"].append(run_metadata)
            
            return metadata
            
        except Exception as e:
            print(f"Error parsing experiment: {e}", file=sys.stderr)
            return None

    def _get_text(self, parent: ET.Element, xpath: str) -> Optional[str]:
        """Safely extract text from XML element."""
        if parent is None:
            return None
        elem = parent.find(xpath)
        return elem.text if elem is not None else None

    def _make_request(self, endpoint: str, params: Dict) -> requests.Response:
        """Make request to NCBI E-utilities with rate limiting."""
        url = f"{self.base_url}/{endpoint}"
        response = self.session.get(url, params=params)
        response.raise_for_status()
        
        # NCBI rate limiting (reduced for faster processing)
        time.sleep(0.1)  # ~10 requests per second
        
        return response

    def extract_metadata(self, query: str, max_results: int = 1000, batch_size: Optional[int] = None) -> Dict:
        """Extract comprehensive metadata for studies matching the query."""
        print(f"Searching for studies matching: {query}", file=sys.stderr)
        bioproject_ids = self.search_studies(query, max_results)
        print(f"Found {len(bioproject_ids)} BioProject IDs", file=sys.stderr)
        
        # Limit processing for testing if batch_size specified
        if batch_size and batch_size < len(bioproject_ids):
            bioproject_ids = bioproject_ids[:batch_size]
            print(f"Processing only first {batch_size} studies for testing", file=sys.stderr)
        
        all_metadata = {
            "query": query,
            "bioproject_count": len(bioproject_ids),
            "studies": []
        }
        
        for i, bioproject_id in enumerate(bioproject_ids):
            try:
                print(f"Processing BioProject {bioproject_id} ({i+1}/{len(bioproject_ids)})", file=sys.stderr)
                
                study_metadata = self.get_study_metadata(bioproject_id)
                experiment_ids = self.get_linked_experiments(bioproject_id)
                experiment_metadata = self.get_experiment_metadata(experiment_ids)
                
                study_data = {
                    "bioproject_metadata": study_metadata,
                    "experiments": experiment_metadata,
                    "experiment_count": len(experiment_metadata),
                    "total_runs": sum(len(exp.get("runs", [])) for exp in experiment_metadata)
                }
                
                all_metadata["studies"].append(study_data)
                
                # Progress checkpoint every 10 studies
                if (i + 1) % 10 == 0:
                    print(f"Completed {i+1}/{len(bioproject_ids)} studies", file=sys.stderr)
                
            except Exception as e:
                print(f"Error processing BioProject {bioproject_id}: {e}", file=sys.stderr)
                continue
                
        return all_metadata


def main():
    parser = argparse.ArgumentParser(
        description="Extract metadata from SRA/ENA for studies matching a query"
    )
    parser.add_argument("query", help="Search query for studies")
    parser.add_argument("-o", "--output", help="Output JSON file", required=True)
    parser.add_argument("--email", help="Email for NCBI API requests")
    parser.add_argument("--max-results", type=int, default=1000,
                        help="Maximum number of studies to process")
    parser.add_argument("--batch-size", type=int, default=None,
                        help="Process only first N studies for testing")
    
    args = parser.parse_args()
    
    extractor = SRAMetadataExtractor(email=args.email)
    metadata = extractor.extract_metadata(args.query, args.max_results, args.batch_size)
    
    with open(args.output, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"Extracted metadata for {metadata['bioproject_count']} studies")
    print(f"Total experiments: {sum(study['experiment_count'] for study in metadata['studies'])}")
    print(f"Total runs: {sum(study['total_runs'] for study in metadata['studies'])}")


if __name__ == "__main__":
    main()