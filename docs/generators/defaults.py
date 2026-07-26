"""
Default pages for pipeline documentation.

These pages are only used when a documentation page does not already
exist.

Existing documentation is NEVER overwritten.
"""

from __future__ import annotations

from .models import Pipeline


README = """# {title}

Describe the purpose of this pipeline.

## Documentation

- Overview
- Input
- Output
- Parameters
- Workflows
- Modules

"""


INDEX = """# {title}

Welcome to the documentation for the **{title}** pipeline.

## Documentation

- [Input](input.md)
- [Output](output.md)
- [Parameters](parameters.md)
- [Workflows](workflows/index.md)
- [Modules](modules/index.md)

"""


INPUT = """# Input

Describe the required input files.

## Required files

| File | Description |
|------|-------------|

"""


OUTPUT = """# Output

Describe the generated output.

## Directory layout

```text
results/