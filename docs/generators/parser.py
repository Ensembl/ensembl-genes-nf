"""
Parser for Ensembl Nextflow pipelines.

The parser builds an in-memory representation of a pipeline that can
later be rendered as Markdown (MkDocs) or MyST/Sphinx.

This module DOES NOT generate documentation.
"""

from __future__ import annotations

import re
from pathlib import Path

from .models import (
    DocumentationPage,
    Directive,
    IOBlock,
    Module,
    Pipeline,
    Workflow,
)

# ---------------------------------------------------------------------
# Nextflow directives recognised by the parser
# ---------------------------------------------------------------------

DIRECTIVES = {
    "label",
    "tag",
    "publishDir",
    "container",
    "conda",
    "cpus",
    "memory",
    "time",
    "cache",
    "storeDir",
    "maxForks",
    "errorStrategy",
}

# ---------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------


def discover_pipelines(root: Path) -> list[Pipeline]:
    """
    Discover every pipeline under pipelines/.

    Expected layout

        pipelines/
            statistics/
            repeats/
            annotation/

    Returns
    -------
    list[Pipeline]
    """

    pipelines: list[Pipeline] = []

    docs_root = root.parent / "docs" / "source" / "pipelines"

    for pipeline_dir in sorted(root.iterdir()):

        if not pipeline_dir.is_dir():
            continue
        print(f"Pipeline: {pipeline_dir.name}")
        pipeline = Pipeline(
            name=pipeline_dir.name,
            root=pipeline_dir,
        )

        #
        # modules
        #

        modules_dir = pipeline_dir / "modules"

        if modules_dir.exists():

            for nf in sorted(modules_dir.glob("*.nf")):

                pipeline.add_module(parse_module(nf))

        #
        # workflow(s)
        #

        workflow = pipeline_dir / "main.nf"

        if workflow.exists():

            pipeline.workflow_file = workflow

            pipeline.add_workflow(parse_workflow(workflow))

        #
        # schema
        #

        schema = pipeline_dir / "nextflow_schema.json"

        if schema.exists():

            pipeline.schema = schema

        #
        # documentation pages
        #

        docs = docs_root / pipeline.name

        pipeline.readme = _page(docs / "README.md")

        pipeline.index = _page(docs / "index.md")

        pipeline.input_doc = _page(docs / "input.md")

        pipeline.output_doc = _page(docs / "output.md")

        pipeline.troubleshooting = _page(docs / "troubleshooting.md")

        pipeline.sort()

        pipelines.append(pipeline)

    return pipelines


def parse_workflow(path: Path) -> Workflow:
    """
    Parse a pipeline workflow (usually main.nf).
    """

    text = path.read_text(
        encoding="utf8",
    )

    workflow = Workflow(
        name=path.stem,
        source=path,
    )

    workflow.description = _clean_comment(_all_comments(text))

    workflow.modules = _discover_workflow_modules(
        text,
    )

    return workflow


def parse_module(path: Path) -> Module:
    """
    Parse one Nextflow module.
    """

    text = path.read_text(
        encoding="utf8",
    )

    module = Module(
        name=path.stem,
        process=_process_name(text),
        source=path,
    )

    #
    # documentation block
    #

    documentation = _parse_documentation(text)

    module.summary = documentation.get(
        "summary",
        "",
    )

    module.description = documentation.get(
        "description",
        "",
    )

    module.tool = documentation.get(
        "tool",
        "",
    )

    module.category = documentation.get(
        "category",
        "",
    )

    module.dependencies = documentation.get(
        "dependencies",
        [],
    )

    module.documented_inputs = documentation.get(
        "inputs",
        [],
    )

    module.documented_outputs = documentation.get(
        "outputs",
        [],
    )

    #
    # Nextflow directives
    #

    lines = text.splitlines()

    _parse_directives(
        module,
        lines,
    )
    module.inputs = _parse_io(
        lines,
        "input:",
    )

    module.outputs = _parse_io(
        lines,
        "output:",
    )

    module.script_summary = _summarise_script(
        lines,
    )
    print("PARSER")
    print(module.process)
    print(module.documented_inputs)
    return module


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def _page(path: Path) -> DocumentationPage | None:
    """
    Return a documentation page if it exists.
    """

    if not path.exists():
        return None

    return DocumentationPage(
        path=path,
        generated=False,
    )


def _process_name(text: str) -> str:

    match = re.search(
        r"^\s*process\s+([A-Za-z0-9_]+)\b",
        text,
        re.MULTILINE,
    )

    if match:
        return match.group(1)

    return "UNKNOWN"


# ---------------------------------------------------------------------
# Documentation block
# ---------------------------------------------------------------------

DOCUMENTATION_TAGS = {
    "tool",
    "category",
    "summary",
    "description",
    "inputs",
    "outputs",
    "dependencies",
}


def _parse_documentation(text: str) -> dict:
    """
    Parse documentation blocks.

    Supports both structured (@summary, @description, ...)
    and legacy comments.
    """

    result: dict = {}

    blocks = re.findall(
        r"/\*(.*?)\*/",
        text,
        re.S,
    )

    if not blocks:
        return result

    #
    # 1. Structured documentation
    #

    for block in blocks:
        if any(
            tag in block
            for tag in (
                "@summary",
                "@description",
                "@inputs",
                "@outputs",
                "@tool",
                "@category",
                "@dependencies",
                "@implementation",
            )
        ):
            return _parse_structured_doc(block)

    #
    # 2. Legacy documentation
    #

    for block in blocks:

        if "Licensed under the Apache License" in block:
            continue

        result["description"] = _clean_comment(block)
        return result

    return result


# ---------------------------------------------------------------------


def _parse_structured_doc(block: str) -> dict:

    result = {}

    current: str | None = None

    values: list[str] = []

    def flush():

        nonlocal current
        nonlocal values

        if current is None:
            return

        value = [v for v in values if v]

        if current in (
            "inputs",
            "outputs",
            "dependencies",
        ):

            result[current] = value

        else:

            result[current] = "\n".join(value)

    for raw in block.splitlines():

        line = raw.strip()

        line = line.lstrip("*").strip()

        if not line:
            continue

        #
        # New section?
        #

        if line.startswith("@"):

            flush()

            tag = line[1:].strip()

            current = tag

            values = []

            continue

        if current:

            values.append(line)

    flush()

    return result


# ---------------------------------------------------------------------


def _clean_comment(comment: str) -> str:

    lines = []

    for raw in comment.splitlines():

        line = raw.strip()

        line = line.lstrip("*").strip()

        if not line:
            continue

        if set(line) <= {"-", "="}:
            continue

        lines.append(line)

    return "\n".join(lines)

# ---------------------------------------------------------------------
# Nextflow directives
# ---------------------------------------------------------------------


def _parse_directives(
    module: Module,
    lines: list[str],
) -> None:
    """
    Parse process directives.

    Example

        label "python"
        publishDir "${params.outdir}"
    """

    for raw in lines:

        line = raw.strip()

        if not line:
            continue

        for directive in DIRECTIVES:

            if not line.startswith(directive):
                continue

            value = line[len(directive) :].strip()

            module.directives.append(
                Directive(
                    name=directive,
                    value=value,
                )
            )

            #
            # Frequently used directives
            #

            if directive == "label":
                module.label = value.strip('"')

            elif directive == "tag":
                module.tag = value.strip('"')

            elif directive == "publishDir":
                module.publish_dir = value

            break


# ---------------------------------------------------------------------
# Input / Output parsing
# ---------------------------------------------------------------------


def _parse_io(
    lines: list[str],
    keyword: str,
) -> list[IOBlock]:
    """
    Parse the Nextflow input/output block.

    Example

        tuple val(meta), path(genome), emit: genome_out

    becomes

        qualifier = tuple

        declaration = val(meta), path(genome)

        emit = genome_out
    """

    block: list[IOBlock] = []

    inside = False

    for raw in lines:

        line = raw.rstrip()

        if line.strip() == keyword:

            inside = True

            continue

        if not inside:
            continue

        stripped = line.strip()

        #
        # End of block?
        #

        if stripped in (
            "input:",
            "output:",
            "script:",
            "exec:",
            "stub:",
            "when:",
        ):
            break

        if not stripped:
            continue

        #
        # qualifier
        #

        match = re.match(
            r"(\w+)\s+(.*)",
            stripped,
        )

        if match is None:

            block.append(
                IOBlock(
                    raw=stripped,
                    qualifier="",
                    declaration=stripped,
                    emit=None,
                )
            )

            continue

        qualifier = match.group(1)

        declaration = match.group(2)

        emit = None

        #
        # emit:
        #

        if ", emit:" in declaration:

            declaration, emit = declaration.rsplit(
                ", emit:",
                1,
            )

            emit = emit.strip()

        block.append(
            IOBlock(
                raw=stripped,
                qualifier=qualifier,
                declaration=declaration.strip(),
                emit=emit,
            )
        )

    return block


# ---------------------------------------------------------------------
# Script summary
# ---------------------------------------------------------------------


def _summarise_script(#pylint: disable=too-many-branches,too-many-statements
    lines: list[str],
) -> list[str]:
    """
    Produce a lightweight summary of what the module does.

    This is heuristic-based.
    """

    summary: list[str] = []

    inside = False

    for raw in lines:

        line = raw.strip()

        if line.startswith("script:"):

            inside = True

            continue

        if not inside:
            continue

        #
        # Wrapper scripts
        #

        if line.startswith("run_"):

            tool = line.split()[0].replace("_", " ").title()

            summary.append(f"Execute {tool}")

        #
        # Python
        #

        elif line.startswith("python "):

            summary.append("Execute Python script")

        #
        # Perl
        #

        elif line.startswith("perl "):

            summary.append("Execute Perl script")

        #
        # Downloads
        #

        elif line.startswith(("wget", "curl")):

            summary.append("Download external resources")

        #
        # File operations
        #

        elif line.startswith("mv "):

            summary.append("Rename output files")

        elif line.startswith("cp "):

            summary.append("Copy output files")

        elif line.startswith("ln -s"):

            summary.append("Create symbolic links")

        #
        # Version report
        #

        elif "versions.yml" in line:

            summary.append("Generate software version report")

    #
    # Remove duplicates preserving order
    #

    seen = set()

    result = []

    for item in summary:

        if item in seen:
            continue

        seen.add(item)

        result.append(item)

    return result


def _all_comments(text: str) -> str:
    """
    Return all /* ... */ comment blocks.
    """

    blocks = re.findall(
        r"/\*(.*?)\*/",
        text,
        flags=re.S,
    )
    print(blocks)
    return "\n\n".join(block.strip() for block in blocks)


def _first_comment(text: str) -> str:
    """
    Return the first /* ... */ block.
    """

    match = re.search(
        r"/\*(.*?)\*/",
        text,
        re.S,
    )

    if match:
        return match.group(1)

    return ""


def _discover_workflow_modules(
    text: str,
) -> list[str]:
    """
    Discover module invocations inside a workflow.

    Example

        FETCH_GENOME(...)
        RUN_REPEATMASKER(...)
    """

    modules = []

    #
    # Matches
    #
    # RUN_REPEATMASKER(
    # FETCH_GENOME(
    #

    pattern = re.compile(
        r"^\s*([A-Z][A-Z0-9_]*)\s*\(",
        re.M,
    )

    for match in pattern.finditer(text):

        module = match.group(1)

        if module not in modules:
            modules.append(module)

    return modules
