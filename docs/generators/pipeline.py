"""
Pipeline documentation generator.

Coordinates the generation of all documentation for a single
Nextflow pipeline.

This module does not render Markdown itself. Rendering is delegated
to specialised generators.
"""

from __future__ import annotations

from pathlib import Path
import logging

from docs.generators.workflow_page import render_pipeline_workflows

from .models import DocumentationPage
from .models import Pipeline
from .module_page import render_pipeline_modules

# from .module_page import render_pipeline_modules
from .parameters import load_schema
from .parameters import render as render_parameters
from .utils import write_file

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------


def generate_pipeline(
    pipeline: Pipeline,
    docs_root: Path,
) -> None:
    """
    Generate the documentation for one pipeline.

    Parameters
    ----------
    pipeline
        Parsed pipeline object.

    docs_root
        Root directory for documentation.
    """

    logger.info(
        "Generating documentation for %s",
        pipeline.name,
    )

    pipeline_dir = docs_root / "pipelines" / pipeline.name

    pipeline_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    #
    # Ensure documentation skeleton exists.
    #

    _ensure_manual_pages(
        pipeline,
        pipeline_dir,
    )

    #
    # Parameter documentation.
    #

    if pipeline.schema:

        logger.info("Generating parameter reference")

        schema = load_schema(
            pipeline.schema,
        )

        write_file(
            pipeline_dir / "parameters.md",
            render_parameters(schema),
        )

    #
    # Module documentation.
    #

    logger.info("Generating module pages")

    render_pipeline_modules(
        pipeline,
        pipeline_dir,
    )

    #
    # Workflow documentation.
    #

    logger.info("Generating workflow pages")

    render_pipeline_workflows(
        pipeline,
        pipeline_dir,
    )


# ---------------------------------------------------------------------
# Manual pages
# ---------------------------------------------------------------------


def _ensure_manual_pages(
    pipeline: Pipeline,
    pipeline_dir: Path,
) -> None:
    """
    Create missing handwritten pages.

    Existing pages are NEVER overwritten.
    """

    pipeline_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    _ensure_page(
        pipeline.readme,
        pipeline_dir / "README.md",
        _readme_template(
            pipeline,
        ),
    )

    _ensure_page(
        pipeline.index,
        pipeline_dir / "index.md",
        _index_template(
            pipeline,
        ),
    )

    _ensure_page(
        pipeline.input_doc,
        pipeline_dir / "input.md",
        "# Input\n",
    )

    _ensure_page(
        pipeline.output_doc,
        pipeline_dir / "output.md",
        "# Output\n",
    )

    _ensure_page(
        pipeline.troubleshooting,
        pipeline_dir / "troubleshooting.md",
        "# Troubleshooting\n",
    )


def _ensure_page(
    page: DocumentationPage | None,
    destination: Path,
    content: str,
) -> None:
    """
    Create a page only if it does not exist.
    """

    if page is not None:
        return

    if destination.exists():
        return

    logger.info(
        "Creating %s",
        destination.name,
    )

    destination.write_text(
        content,
        encoding="utf8",
    )


# ---------------------------------------------------------------------
# Templates
# ---------------------------------------------------------------------


def _readme_template(
    pipeline: Pipeline,
) -> str:

    return f"""# {pipeline.title}

Welcome to the documentation for the
**{pipeline.title}** pipeline.

## Documentation

```{{toctree}}
:maxdepth: 1

input
output
parameters
modules/index
workflows/index
troubleshooting
```
"""


def _index_template(
    pipeline: Pipeline,
) -> str:

    return f"""# {pipeline.title}

Welcome to the documentation for the
**{pipeline.title}** pipeline.

```{{toctree}}
:maxdepth: 1

README
input
output
parameters
modules/index
workflows/index
troubleshooting
```
"""
