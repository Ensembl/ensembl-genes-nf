#!/usr/bin/env python3
"""
Generate documentation for all Ensembl Nextflow pipelines.

The generator performs the following steps:

1. Discover pipelines
2. Parse Nextflow modules
3. Generate module documentation
4. Generate workflow documentation
5. Generate parameter reference
6. Generate navigation pages

Only generated documentation is overwritten.
Handwritten documentation is preserved.
"""

from __future__ import annotations

import logging
from pathlib import Path

#from generators.navigation import generate_navigation
#from generators.parser import discover_pipelines
#from generators.pipeline import generate_pipeline
from docs.generators.parser import discover_pipelines
from docs.generators.pipeline import generate_pipeline
from docs.generators.navigation import (
    generate_navigation,
    generate_root_index,
)
# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------
# Repository layout
# ---------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]

PIPELINES_DIR = REPO_ROOT / "pipelines"

DOCS_DIR = REPO_ROOT / "docs" / "source"


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main() -> None:
    """
    Generate the documentation for every pipeline.
    """

    logger.info("Discovering pipelines...")

    pipelines = discover_pipelines(
        PIPELINES_DIR,
    )

    logger.info(
        "Found %d pipeline(s).",
        len(pipelines),
    )

    for pipeline in pipelines:

        logger.info(
            "Generating documentation for '%s'...",
            pipeline.name,
        )

        generate_pipeline(
            pipeline=pipeline,
            docs_root=DOCS_DIR,
        )

    logger.info(
        "Generating documentation index..."
    )

    generate_navigation(
        docs_root=DOCS_DIR,
        pipelines=pipelines,
    )
    generate_root_index(
        docs_root=DOCS_DIR,
    )

    logger.info("Documentation completed successfully.")


# ---------------------------------------------------------------------

if __name__ == "__main__":
    main()