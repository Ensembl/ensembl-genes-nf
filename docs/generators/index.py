"""
Generation of documentation indexes.

This module creates the landing pages for the generated documentation.
"""

from __future__ import annotations

from pathlib import Path

from .markdown import write
from .models import Pipeline


def render_root(pipelines: list[Pipeline]) -> str:
    """
    Generate the root generated/index.md page.
    """

    out = [
        "# Generated Documentation",
        "",
        "These pages are generated automatically from the Nextflow pipelines.",
        "",
    ]

    for pipeline in sorted(pipelines, key=lambda p: p.name):

        out.extend(
            [
                f"## {pipeline.title}",
                "",
                f"- [{pipeline.title} modules]({pipeline.name}/index.md)",
                f"- [{pipeline.title} parameters](../generated/{pipeline.name}-parameters.md)",
                "",
            ]
        )

    return "\n".join(out)


def render_pipeline(pipeline: Pipeline) -> str:
    """
    Landing page for one pipeline.
    """

    out = [
        f"# {pipeline.title}",
        "",
        f"Automatically generated documentation for the **{pipeline.title}** pipeline.",
        "",
        "## Modules",
        "",
        "| Process | Module | Description |",
        "|---------|--------|-------------|",
    ]

    for module in pipeline.modules:

        description = module.description.splitlines()[0] if module.description else ""

        out.append(
            "| "
            f"`{module.process}` | "
            f"[{module.slug}](modules/{module.slug}.md)"
            f" | {description} |"
        )

    out.append("")

    return "\n".join(out)


def write_indexes(
    generated_dir: Path,
    pipelines: list[Pipeline],
) -> None:
    """
    Write all index pages.
    """

    #
    # root index
    #

    write(
        generated_dir / "index.md",
        render_root(pipelines),
    )

    #
    # pipeline indexes
    #

    for pipeline in pipelines:

        write(
            generated_dir / pipeline.name / "index.md",
            render_pipeline(pipeline),
        )