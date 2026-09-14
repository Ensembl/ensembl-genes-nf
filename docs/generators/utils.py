"""
Common utilities for the documentation generator.
"""

from __future__ import annotations

import re
from pathlib import Path


# ---------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------


def ensure_dir(path: Path) -> None:
    """Create a directory if it does not already exist."""
    path.mkdir(parents=True, exist_ok=True)


def write_file(path: Path, text: str) -> None:
    """Write a UTF-8 text file."""
    ensure_dir(path.parent)
    path.write_text(text.rstrip() + "\n", encoding="utf8")


# ---------------------------------------------------------------------
# Strings
# ---------------------------------------------------------------------


def slugify(name: str) -> str:
    """
    Convert a module name into a documentation slug.

    Example:
        run_repeatmasker -> run-repeatmasker
    """
    return name.replace("_", "-")


def title_case(name: str) -> str:
    """
    Convert an identifier into a human-readable title.

    Example:
        repeat_annotation -> Repeat Annotation
    """
    return name.replace("_", " ").title()


def clean_comment(text: str) -> str:
    """
    Remove comment markers and empty lines from a Nextflow block comment.
    """

    lines = []

    for line in text.splitlines():

        line = line.strip()

        line = line.lstrip("*").strip()

        if not line:
            continue

        if line.startswith("="):
            continue

        if line.startswith("---"):
            continue

        if line.startswith("Licensed"):
            continue

        if line.startswith("See the NOTICE"):
            continue

        if line.startswith("http"):
            continue

        lines.append(line)

    return "\n".join(lines)


# ---------------------------------------------------------------------
# Regex helpers
# ---------------------------------------------------------------------


PROCESS_RE = re.compile(r"process\s+([A-Za-z0-9_]+)")
COMMENT_RE = re.compile(r"/\*(.*?)\*/", re.S)


def first_comment(text: str) -> str:
    """
    Return the first /* ... */ block.
    """

    match = COMMENT_RE.search(text)

    if match is None:
        return ""

    return clean_comment(match.group(1))


def first_process(text: str) -> str:
    """
    Return the first process name.
    """

    match = PROCESS_RE.search(text)

    if match is None:
        return "UNKNOWN"

    return match.group(1)


# ---------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------


def md_heading(title: str, level: int = 1) -> str:#pylint: disable=missing-function-docstring
    return "#" * level + f" {title}"


def md_table(headers: list[str]) -> list[str]:
    """
    Create a markdown table header.

    Example:

    md_table(["Name","Value"])
    """

    return [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("-" * (len(h) + 2) for h in headers) + "|",
    ]


def md_code(lines: list[str], language: str = "text") -> list[str]:
    """
    Create a fenced code block.
    """

    return [f"```{language}"] + lines + ["```"]
