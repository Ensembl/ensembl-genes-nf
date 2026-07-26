"""
Markdown renderer for generated documentation.
"""

from __future__ import annotations

from importlib.resources import path
from pathlib import Path

from .models import Module
from .models import Pipeline


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------


def write(path: Path, text: str) -> None:
    """Write file creating parents if necessary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf8")


def heading(title: str, level: int = 1) -> str:
    return "#" * level + f" {title}"


# ---------------------------------------------------------------------
# Module pages
# ---------------------------------------------------------------------

from pathlib import Path

from .models import Pipeline
from .utils import write_file


def render_pipeline_modules(
    pipeline: Pipeline,
    pipeline_dir: Path,
) -> None:
    """
    Generate all module pages for one pipeline.
    """

    if not pipeline.modules:
        return

    modules_dir = pipeline_dir / "modules"

    modules_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    #
    # Individual module pages
    #

    for module in pipeline.modules:

        write_file(
            modules_dir / f"{module.slug}.md",
            render_module(module),
        )

    #
    # modules/index.md
    #

    lines = [
        "# Modules",
        "",
        "Documentation for the modules used by this pipeline.",
        "",
        "```{toctree}",
        ":maxdepth: 1",
        "",
    ]

    for module in sorted(
        pipeline.modules,
        key=lambda m: m.slug,
    ):
        lines.append(module.slug)

    lines.extend([
        "```",
        "",
    ])

    write_file(
        modules_dir / "index.md",
        "\n".join(lines),
    )
    
def render_module(module: Module) -> str:
    """
    Render one module page.
    """

    lines: list[str] = []

    lines.append(heading(module.process))
    lines.append("")

    if module.description:
        lines.append(module.description)
        lines.append("")

    #
    # Process
    #

    lines.append(heading("Process Details", 2))
    lines.append("")

    lines.append("| Property | Value |")
    lines.append("|----------|-------|")

    lines.append(f"| Process | `{module.process}` |")

    if module.label:
        lines.append(f"| Label | `{module.label}` |")

    if module.tag:
        lines.append(f"| Tag | `{module.tag}` |")

    if module.publish_dir:
        lines.append(f"| Publish directory | `{module.publish_dir}` |")

    for directive in module.directives:

        if directive.name in {"label", "tag", "publishDir"}:
            continue

        lines.append(
            f"| {directive.name} | `{directive.value}` |"
        )

    lines.append("")

    #
    # Inputs
    #

    lines.append(heading("Inputs", 2))
    lines.append("")

    if module.documented_inputs:

        lines.append("| Name | Description |")
        lines.append("| ---- | ----------- |")

        for entry in module.documented_inputs:
            name, _, description = entry.partition(" ")
            lines.append(f"| `{name}` | {description.strip()} |")

        lines.append("")

    if module.inputs:

        lines.append("### Nextflow interface")
        lines.append("")
        lines.append("```nextflow")

        for item in module.inputs:
            lines.append(item.raw)

        lines.append("```")

    if not module.documented_inputs and not module.inputs:

        lines.append("*No inputs documented.*")

    lines.append("")

    #
    # Outputs
    #

    lines.append(heading("Outputs", 2))
    lines.append("")

    if module.documented_outputs:

        lines.append("| Name | Description |")
        lines.append("| ---- | ----------- |")

        for entry in module.documented_outputs:
            name, _, description = entry.partition(" ")
            lines.append(f"| `{name}` | {description.strip()} |")

        lines.append("")

    if module.outputs:

        lines.append("### Nextflow interface")
        lines.append("")
        lines.append("```nextflow")

        for item in module.outputs:
            lines.append(item.raw)

        lines.append("```")

    if not module.documented_outputs and not module.outputs:

        lines.append("*No outputs documented.*")

    lines.append("")

    #
    # Implementation
    #

    if module.script_summary:

        lines.append(heading("Implementation Summary", 2))
        lines.append("")

        for step in module.script_summary:
            lines.append(f"- {step}")

        lines.append("")

    #
    # Source
    #

    lines.append(heading("Source", 2))
    lines.append("")

    lines.append(
        f"`{module.source}`"
    )

    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------
# Pipeline index
# ---------------------------------------------------------------------


def render_pipeline_index(pipeline: Pipeline) -> str:

    lines = []

    lines.append(heading(f"{pipeline.title} Modules"))
    lines.append("")

    lines.append(
        "This page is generated automatically from the Nextflow modules."
    )

    lines.append("")

    lines.append("| Process | Description |")

    lines.append("|---------|-------------|")

    for module in pipeline.modules:

        description = module.description.split("\n")[0]

        lines.append(
            "| "
            f"[`{module.process}`](modules/{module.slug}.md)"
            f" | {description} |"
        )

    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------
# Global index
# ---------------------------------------------------------------------


def render_global_index(
    pipelines: list[Pipeline],
) -> str:

    lines = []

    lines.append(heading("Generated Documentation"))

    lines.append("")

    lines.append(
        "The pages in this directory are generated automatically."
    )

    lines.append("")

    for pipeline in pipelines:

        lines.append(heading(pipeline.title, 2))
        lines.append("")

        lines.append(
            f"- [{pipeline.title} modules]({pipeline.name}/index.md)"
        )

        lines.append(
            f"- [{pipeline.title} parameters](../generated/{pipeline.name}-parameters.md)"
        )

        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------
# Write helpers
# ---------------------------------------------------------------------


def write_pipeline(
    root: Path,
    pipeline: Pipeline,
) -> None:

    pipeline_dir = root / pipeline.name

    write(
        pipeline_dir / "index.md",
        render_pipeline_index(pipeline),
    )

    modules_dir = pipeline_dir / "modules"

    for module in pipeline.modules:

        write(
            modules_dir / f"{module.slug}.md",
            render_module(module),
        )


def write_global(
    root: Path,
    pipelines: list[Pipeline],
) -> None:

    write(
        root / "index.md",
        render_global_index(pipelines),
    )