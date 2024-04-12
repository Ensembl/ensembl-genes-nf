# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
import sys
import datetime

sys.path.insert(0, os.path.abspath("../../src/python/ensembl/genes/metadata"))

print(sys.executable)

# -- Project information -----------------------------------------------------
project = 'ensembl-genes-metadata'
author = 'ensembl@dev.org'
copyright_owner = "EMBL-European Bioinformatics Institute"
copyright_dates = "[2016-%d]" % datetime.datetime.now().year
copyright = copyright_dates + " " + copyright_owner
html_baseurl = 'https://ensembl.github.io/ensembl-genes-metadata/'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = "0.1"
# The full version, including alpha/beta/rc tags.
release = "0.1"

copyright_owner = "EMBL-European Bioinformatics Institute"
copyright_dates = "[2016-%d]" % datetime.datetime.now().year
copyright = copyright_dates + " " + copyright_owner
html_baseurl = 'https://ensembl.github.io/ensembl-genes-metadata/'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    
]

# Defining autodoc functionality
autodoc_default_options = {
    "member-order": "alphabetical",
    "undoc-members": False,
    "exclude-members": "__weakref__",
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

source_suffix = {".rst": "restructuredtext"}

# The master toctree document.
master_doc = "index"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# A list of ignored prefixes for module index sorting.
modindex_common_prefix = ["ensembl."]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "agogo"
html_theme_options = {
    "bodyfont": "Garamond, Arial, serif",
    "headerfont": "Arial, Helvetica, serif",
    "headerlinkcolor": "#33d6ff",
    "pagewidth": "70em",
    "documentwidth": "50em",
    "rightsidebar": True,
    "bgcolor": "#009999",
    "headerbg": "#009999",
    "footerbg": "#e6fff9",
    "linkcolor": "green",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Ensembldoc'
# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
#man_pages = [(master_doc, "ensembl-genes-metadata", "Ensembl Anno Base Library Documentation", [author], 1)]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "ensembl-genes-metadata.tex", "Ensembl-genes-metadata Documentation", [author], "manual"),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "ensembl-genes-metadata.tex", "Ensembl-genes-metadata Documentation", [author], "manual"),]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "ensembl-genes-metadata",
        "Ensembl-genes-metadata Documentation",
        author,
        "ensembl-genes-metadata",
        "Ensembl-genes-metadata Documentation.",
        "Miscellaneous",
    ),
]