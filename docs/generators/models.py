"""
Data models used by the documentation generator.

These dataclasses represent the documentation objects extracted from
Nextflow pipelines and JSON schemas.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------
# Nextflow process information
# ---------------------------------------------------------------------


@dataclass(slots=True)
class Directive:
    """A Nextflow process directive."""

    name: str
    value: str


@dataclass(slots=True)
class IOBlock:
    """
    Parsed Nextflow input/output declaration.

    Example
    -------
    tuple val(meta), path(genome), emit: genome_out
    """

    raw: str

    qualifier: str

    declaration: str

    emit: str | None = None


# ---------------------------------------------------------------------
# Documentation pages
# ---------------------------------------------------------------------


@dataclass(slots=True)
class DocumentationPage:
    """
    Represents one documentation page.

    generated=False -> handwritten page.
    generated=True  -> generated automatically.
    """

    path: Path

    title: str = ""

    generated: bool = False


# ---------------------------------------------------------------------
# Nextflow module
# ---------------------------------------------------------------------


@dataclass(slots=True)
class Module:
    """Representation of one Nextflow module."""

    #
    # identity
    #

    name: str
    process: str
    source: Path

    pipeline: str = ""

    #
    # documentation
    #

    summary: str = ""

    description: str = ""

    category: str = ""

    tool: str = ""

    dependencies: list[str] = field(default_factory=list)

    documented_inputs: list[str] = field(default_factory=list)

    documented_outputs: list[str] = field(default_factory=list)

    #
    # nextflow directives
    #

    label: str | None = None

    tag: str | None = None

    publish_dir: str | None = None

    directives: list[Directive] = field(default_factory=list)

    #
    # nextflow IO
    #

    inputs: list[IOBlock] = field(default_factory=list)

    outputs: list[IOBlock] = field(default_factory=list)

    #
    # implementation summary
    #

    script_summary: list[str] = field(default_factory=list)

    @property
    def slug(self) -> str:
        return self.name.replace("_", "-")

    @property
    def title(self) -> str:
        return self.process


# ---------------------------------------------------------------------
# Workflow
# ---------------------------------------------------------------------


@dataclass(slots=True)
class Workflow:
    """
    Representation of a workflow.

    Usually corresponds to main.nf.
    """

    name: str

    source: Path

    description: str = ""

    modules: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------


@dataclass(slots=True)
class Pipeline:
    """Representation of a complete pipeline."""

    name: str

    root: Path

    workflow_file: Path | None = None

    schema: Path | None = None

    modules: list[Module] = field(default_factory=list)

    workflows: list[Workflow] = field(default_factory=list)

    #
    # handwritten documentation
    #

    readme: DocumentationPage | None = None

    index: DocumentationPage | None = None

    input_doc: DocumentationPage | None = None

    output_doc: DocumentationPage | None = None

    troubleshooting: DocumentationPage | None = None

    @property
    def title(self) -> str:
        return self.name.replace("_", " ").title()

    @property
    def slug(self) -> str:
        return self.name.replace("_", "-")

    def add_module(self, module: Module) -> None:
        module.pipeline = self.name
        self.modules.append(module)

    def add_workflow(self, workflow: Workflow) -> None:
        self.workflows.append(workflow)

    def sort(self) -> None:
        self.modules.sort(key=lambda m: m.process)
        self.workflows.sort(key=lambda w: w.name)


# ---------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------


@dataclass(slots=True)
class Parameter:
    """Single parameter from nextflow_schema.json."""

    name: str

    type: str

    description: str

    default: str

    required: bool


@dataclass(slots=True)
class ParameterGroup:
    """Logical parameter group."""

    title: str

    description: str

    parameters: list[Parameter] = field(default_factory=list)


@dataclass(slots=True)
class Schema:
    """Representation of one nextflow_schema.json."""

    pipeline: str

    groups: list[ParameterGroup] = field(default_factory=list)