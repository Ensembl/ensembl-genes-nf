# pylint: disable=too-many-instance-attributes,missing-function-docstring,too-many-public-methods
"""
Generate parameter reference documentation from Nextflow schemas.
"""

from __future__ import annotations

import json
from pathlib import Path

from .models import Parameter
from .models import ParameterGroup
from .models import Schema


# ---------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------


def load_schema(path: Path) -> Schema:

    with path.open(encoding="utf8") as fh:
        raw = json.load(fh)

    schema = Schema(
        pipeline=path.parent.name,
    )

    defs = raw.get("$defs", {})

    #
    # grouped schema
    #

    for item in raw.get("allOf", []):

        ref = item.get("$ref")

        if not ref:
            continue

        name = ref.split("/")[-1]

        definition = defs[name]

        group = ParameterGroup(
            title=definition.get("title", name),
            description=definition.get("description", ""),
        )

        required = set(definition.get("required", []))

        for pname, prop in definition.get("properties", {}).items():

            group.parameters.append(
                Parameter(
                    name=pname,
                    type=_parameter_type(prop),
                    description=prop.get("description", ""),
                    default=str(prop.get("default", "")),
                    required=pname in required,
                )
            )

        schema.groups.append(group)

    #
    # flat schema
    #

    if not schema.groups:

        group = ParameterGroup(
            title="Parameters",
            description=raw.get("description", ""),
        )

        required = set(raw.get("required", []))

        for pname, prop in raw.get("properties", {}).items():

            group.parameters.append(
                Parameter(
                    name=pname,
                    type=_parameter_type(prop),
                    description=prop.get("description", ""),
                    default=str(prop.get("default", "")),
                    required=pname in required,
                )
            )

        schema.groups.append(group)

    return schema


# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------


def _parameter_type(prop):

    ptype = prop.get("type", "")

    if isinstance(ptype, list):
        ptype = ", ".join(ptype)

    if prop.get("enum"):

        values = ", ".join(map(str, prop["enum"]))

        ptype += f" ({values})"

    return ptype


# ---------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------


def render(schema: Schema) -> str:

    out = []

    out.append(f"# {schema.pipeline.title()} Parameters")
    out.append("")
    out.append("Automatically generated from the pipeline " "`nextflow_schema.json`.")#pylint: disable=implicit-str-concat
    out.append("")

    for group in schema.groups:

        out.append(f"## {group.title}")
        out.append("")

        if group.description:

            out.append(group.description)
            out.append("")

        out.append("| Parameter | Type | Default | Required | Description |")

        out.append("|-----------|------|---------|----------|-------------|")

        for parameter in group.parameters:

            out.append(
                "| "
                f"`{parameter.name}` | "
                f"{parameter.type} | "
                f"{parameter.default} | "
                f"{'yes' if parameter.required else 'no'} | "
                f"{parameter.description} |"
            )

        out.append("")

    return "\n".join(out)
