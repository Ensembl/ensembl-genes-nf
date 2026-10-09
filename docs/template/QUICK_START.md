# Quick start

The example pipeline is the executable starting point for a new pipeline.
It uses Nextflow's strict/v2 syntax. Nextflow 26.04 and newer use the v2
parser by default; Nextflow 25.10.2 needs `NXF_SYNTAX_PARSER=v2`.

## Run the example

```bash
cd pipelines/example
nextflow lint -o concise .
nextflow run main.nf -stub-run -profile test
```

The test profile uses the deterministic input in `assets/input.txt`, runs two
sample records through two tools, combines their outputs, and counts lines. The
outputs are written below `results/test`.

For the smallest example:

```bash
nextflow run workflows/simple_workflow.nf -stub-run --outdir results/simple
```

## Copy the structure

Use this layout for a real pipeline:

```text
pipelines/my-pipeline/
├── main.nf                 # validate parameters and call the entry workflow
├── nextflow.config         # plugins, defaults, and test/profile settings
├── nextflow_schema.json    # parameter contract
├── assets/                 # small deterministic test inputs and references
├── workflows/              # orchestration and conditional branching
├── subworkflows/           # reusable chains of modules
└── modules/                # one process/tool per file
```

Start a module with `modules/module_example.nf`. Use
`modules/reference_example.nf` when you need dynamic arguments, optional
outputs, and retry configuration.

As a default, keep named parameters in `main.nf`. Validate and normalize them
there, then pass explicit channels and values into workflows and subworkflows.
This can keep reusable components independent of one pipeline's parameter names
and schema. It is a design choice rather than a restriction: direct parameter
access may be fine for genuinely pipeline-wide settings.
See [workflows.md](workflows.md#keep-parameter-ownership-at-the-entrypoint)
for the rationale and a concrete example.

## Development loop

1. Add or change the schema and defaults.
2. Update the README and test profile.
3. Keep orchestration in workflows and tool commands in modules.
4. Make the stub create every declared output.
5. Run strict lint and the smallest stub workflow.

```bash
nextflow lint -o concise .
nextflow run main.nf -stub-run -profile test
```

See [NEXTFLOW_REQUIREMENTS.md](../NEXTFLOW_REQUIREMENTS.md) for the repository
contract and [PATTERNS.md](PATTERNS.md) for reusable dataflow patterns.
