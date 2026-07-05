#!/usr/bin/env python3
"""Render small Markdown parameter references from Nextflow JSON schemas."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = REPO_ROOT / "docs" / "source"
GENERATED_DIR = SOURCE_DIR / "generated"

SCHEMAS = {
    "example": REPO_ROOT / "pipelines" / "example" / "nextflow_schema.json",
    "repeat": REPO_ROOT / "pipelines" / "repeat" / "nextflow_schema.json",
    "riboseq": REPO_ROOT / "pipelines" / "riboseq" / "nextflow_schema.json",
    "translon-consensus": REPO_ROOT
    / "pipelines"
    / "translon-consensus"
    / "nextflow_schema.json",
}


def load_schema(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def escape_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return json.dumps(value)
    if isinstance(value, (dict, list)):
        text = json.dumps(value, sort_keys=True)
    else:
        text = str(value)
    return text.replace("|", "\\|").replace("\n", " ")


def property_type(details: dict[str, Any]) -> str:
    value = details.get("type", "")
    if isinstance(value, list):
        value = ", ".join(value)
    enum = details.get("enum")
    if enum:
        value = f"{value}; one of: {', '.join(map(str, enum))}"
    return str(value)


def schema_groups(schema: dict[str, Any]) -> list[tuple[str, str, dict[str, Any], set[str]]]:
    groups: list[tuple[str, str, dict[str, Any], set[str]]] = []

    for item in schema.get("allOf", []):
        ref = item.get("$ref", "")
        if not ref.startswith("#/$defs/"):
            continue
        name = ref.rsplit("/", 1)[-1]
        definition = schema.get("$defs", {}).get(name, {})
        groups.append(
            (
                definition.get("title", name.replace("_", " ").title()),
                definition.get("description", ""),
                definition.get("properties", {}),
                set(definition.get("required", [])),
            )
        )

    if not groups:
        groups.append(
            (
                "Parameters",
                schema.get("description", ""),
                schema.get("properties", {}),
                set(schema.get("required", [])),
            )
        )

    return groups


def render_schema(slug: str, schema_path: Path) -> str:
    schema = load_schema(schema_path)
    title = schema.get("title", f"{slug} parameters")
    lines = [
        f"# {title}",
        "",
        f"Generated from `{schema_path.relative_to(REPO_ROOT)}`.",
        "",
    ]

    for group_title, group_description, properties, required in schema_groups(schema):
        lines.extend([f"## {group_title}", ""])
        if group_description:
            lines.extend([group_description, ""])

        lines.extend(
            [
                "| Parameter | Type | Default | Required | Description |",
                "| --- | --- | --- | --- | --- |",
            ]
        )

        for name, details in properties.items():
            lines.append(
                "| `{}` | {} | {} | {} | {} |".format(
                    escape_cell(name),
                    escape_cell(property_type(details)),
                    escape_cell(details.get("default", "")),
                    "yes" if name in required else "no",
                    escape_cell(details.get("description", "")),
                )
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def render_index(slugs: list[str]) -> str:
    lines = [
        "# Generated Parameter Reference",
        "",
        "These pages are generated from the pipeline `nextflow_schema.json` files.",
        "",
    ]
    lines.extend(f"- [{slug} parameters]({slug}-parameters.md)" for slug in slugs)
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    for slug, schema_path in SCHEMAS.items():
        output_path = GENERATED_DIR / f"{slug}-parameters.md"
        output_path.write_text(render_schema(slug, schema_path), encoding="utf-8")
    (GENERATED_DIR / "index.md").write_text(render_index(list(SCHEMAS)), encoding="utf-8")


if __name__ == "__main__":
    main()
