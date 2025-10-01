#!/usr/bin/env python3

import argparse
import json
import sys
import time
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET
import re

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class PMCContentExtractor:
    def __init__(self, email: Optional[str] = None, max_retries: int = 3):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.pmc_base_url = "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi"
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

    def search_pmc_articles(self, query: str, max_results: int = 100) -> List[str]:
        """Search PMC for articles matching query and return PMC IDs."""
        params = {
            "db": "pmc",
            "term": f"{query}[Title/Abstract]",
            "retmax": max_results,
            "retmode": "xml"
        }
        if self.email:
            params["email"] = self.email
            
        response = self._make_request("esearch.fcgi", params)
        root = ET.fromstring(response.content)
        
        pmc_ids = []
        for id_elem in root.findall(".//Id"):
            pmc_ids.append(id_elem.text)
            
        return pmc_ids

    def get_article_metadata(self, pmc_id: str) -> Dict:
        """Fetch article metadata from PMC."""
        params = {
            "db": "pmc",
            "id": pmc_id,
            "retmode": "xml"
        }
        if self.email:
            params["email"] = self.email
            
        try:
            response = self._make_request("efetch.fcgi", params)
            root = ET.fromstring(response.content)
            
            article = root.find(".//article")
            if article is None:
                return {}
                
            metadata = {
                "pmc_id": f"PMC{pmc_id}",
                "title": self._get_text(article, ".//article-title"),
                "abstract": self._get_text(article, ".//abstract"),
                "journal": self._get_text(article, ".//journal-title"),
                "publication_date": self._extract_date(article),
                "authors": self._extract_authors(article),
                "keywords": self._extract_keywords(article),
                "methods_section": None,
                "full_text_available": True
            }
            
            return metadata
            
        except Exception as e:
            print(f"Error fetching PMC article {pmc_id}: {e}", file=sys.stderr)
            return {"pmc_id": f"PMC{pmc_id}", "full_text_available": False}

    def get_full_text_content(self, pmc_id: str) -> Dict:
        """Extract full text content from PMC article."""
        try:
            # Use PMC OAI service to get full text XML
            params = {
                "verb": "GetRecord",
                "identifier": f"oai:pubmedcentral.nih.gov:{pmc_id}",
                "metadataPrefix": "pmc"
            }
            
            response = self.session.get(self.pmc_base_url, params=params)
            response.raise_for_status()
            
            root = ET.fromstring(response.content)
            article = root.find(".//{http://www.ncbi.nlm.nih.gov/pmc/articleset/}article")
            
            if article is None:
                return {"methods": None, "results": None, "error": "No full text available"}
            
            content = {
                "methods": self._extract_methods_section(article),
                "results": self._extract_results_section(article),
                "materials": self._extract_materials_section(article),
                "experimental_procedures": self._extract_experimental_procedures(article)
            }
            
            return content
            
        except Exception as e:
            print(f"Error fetching full text for PMC{pmc_id}: {e}", file=sys.stderr)
            return {"methods": None, "results": None, "error": str(e)}

    def _extract_methods_section(self, article: ET.Element) -> Optional[str]:
        """Extract methods/materials section from article XML."""
        methods_patterns = [
            ".//sec[@sec-type='methods']",
            ".//sec[@sec-type='materials-methods']", 
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'method')]]",
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'material')]]"
        ]
        
        for pattern in methods_patterns:
            methods_sec = article.find(pattern)
            if methods_sec is not None:
                return self._extract_section_text(methods_sec)
                
        return None

    def _extract_results_section(self, article: ET.Element) -> Optional[str]:
        """Extract results section from article XML."""
        results_patterns = [
            ".//sec[@sec-type='results']",
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'result')]]"
        ]
        
        for pattern in results_patterns:
            results_sec = article.find(pattern)
            if results_sec is not None:
                return self._extract_section_text(results_sec)
                
        return None

    def _extract_materials_section(self, article: ET.Element) -> Optional[str]:
        """Extract materials and methods subsections."""
        materials_patterns = [
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'rna')]]",
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'sequencing')]]",
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'library')]]"
        ]
        
        sections = []
        for pattern in materials_patterns:
            for sec in article.findall(pattern):
                sections.append(self._extract_section_text(sec))
                
        return "\n\n".join(filter(None, sections)) if sections else None

    def _extract_experimental_procedures(self, article: ET.Element) -> Optional[str]:
        """Extract experimental procedures from various section types."""
        proc_patterns = [
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'procedure')]]",
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'protocol')]]",
            ".//sec[title[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'experimental')]]"
        ]
        
        sections = []
        for pattern in proc_patterns:
            for sec in article.findall(pattern):
                sections.append(self._extract_section_text(sec))
                
        return "\n\n".join(filter(None, sections)) if sections else None

    def _extract_section_text(self, section: ET.Element) -> str:
        """Extract clean text from a section element."""
        text_parts = []
        
        # Get section title
        title = section.find(".//title")
        if title is not None and title.text:
            text_parts.append(title.text.strip())
        
        # Get all paragraphs
        for p in section.findall(".//p"):
            if p.text:
                text_parts.append(p.text.strip())
                
        return "\n".join(text_parts)

    def _extract_date(self, article: ET.Element) -> Optional[str]:
        """Extract publication date."""
        date_elem = article.find(".//pub-date[@pub-type='ppub']")
        if date_elem is None:
            date_elem = article.find(".//pub-date[@pub-type='epub']")
        if date_elem is None:
            date_elem = article.find(".//pub-date")
            
        if date_elem is not None:
            year = self._get_text(date_elem, ".//year")
            month = self._get_text(date_elem, ".//month")
            day = self._get_text(date_elem, ".//day")
            
            if year:
                date_str = year
                if month:
                    date_str += f"-{month.zfill(2)}"
                    if day:
                        date_str += f"-{day.zfill(2)}"
                return date_str
                
        return None

    def _extract_authors(self, article: ET.Element) -> List[str]:
        """Extract author names."""
        authors = []
        for author in article.findall(".//contrib[@contrib-type='author']"):
            given_names = self._get_text(author, ".//given-names")
            surname = self._get_text(author, ".//surname")
            
            if surname:
                name = surname
                if given_names:
                    name = f"{given_names} {surname}"
                authors.append(name)
                
        return authors

    def _extract_keywords(self, article: ET.Element) -> List[str]:
        """Extract article keywords."""
        keywords = []
        for kwd in article.findall(".//kwd"):
            if kwd.text:
                keywords.append(kwd.text.strip())
                
        return keywords

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
        
        # NCBI rate limiting
        time.sleep(0.34)  # ~3 requests per second
        
        return response

    def extract_articles_metadata(self, query: str, max_results: int = 100, 
                                 include_full_text: bool = False) -> Dict:
        """Extract comprehensive metadata from PMC articles."""
        print(f"Searching PMC for articles matching: {query}", file=sys.stderr)
        pmc_ids = self.search_pmc_articles(query, max_results)
        print(f"Found {len(pmc_ids)} PMC articles", file=sys.stderr)
        
        all_metadata = {
            "query": query,
            "article_count": len(pmc_ids),
            "articles": []
        }
        
        for pmc_id in pmc_ids:
            try:
                print(f"Processing PMC{pmc_id}", file=sys.stderr)
                
                article_metadata = self.get_article_metadata(pmc_id)
                
                if include_full_text and article_metadata.get("full_text_available"):
                    full_text = self.get_full_text_content(pmc_id)
                    article_metadata.update(full_text)
                
                all_metadata["articles"].append(article_metadata)
                
            except Exception as e:
                print(f"Error processing PMC{pmc_id}: {e}", file=sys.stderr)
                continue
                
        return all_metadata


def main():
    parser = argparse.ArgumentParser(
        description="Extract content from PMC articles matching a query"
    )
    parser.add_argument("query", help="Search query for PMC articles")
    parser.add_argument("-o", "--output", help="Output JSON file", required=True)
    parser.add_argument("--email", help="Email for NCBI API requests")
    parser.add_argument("--max-results", type=int, default=100,
                        help="Maximum number of articles to process")
    parser.add_argument("--full-text", action="store_true",
                        help="Extract full text content (methods, results sections)")
    
    args = parser.parse_args()
    
    extractor = PMCContentExtractor(email=args.email)
    metadata = extractor.extract_articles_metadata(
        args.query, 
        args.max_results, 
        args.full_text
    )
    
    with open(args.output, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"Extracted metadata for {metadata['article_count']} PMC articles")


if __name__ == "__main__":
    main()