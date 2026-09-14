"""Sphinx configuration for the Ensembl Genes Nextflow documentation."""

from __future__ import annotations

from datetime import date
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

project = "Ensembl Genes Nextflow Pipelines"
author = "Ensembl Genebuild"
copyright = f"{date.today().year}, Ensembl"

extensions = [
    "myst_parser",
    "sphinx.ext.extlinks",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]
myst_heading_anchors = 3

html_theme = "sphinx_rtd_theme"
html_title = "Ensembl Genes Nextflow Pipelines"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
}

extlinks = {
    "repo": ("https://github.com/Ensembl/ensembl-genes-nf/blob/main/%s", "%s"),
}

linkcheck_ignore = [
    r"https://github.com/Ensembl/ensembl-genes-nf/.*",
]
