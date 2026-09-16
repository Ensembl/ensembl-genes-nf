from pathlib import Path
import os


PIPELINE = Path(__file__).parents[1]
ROOT_CONFIG = PIPELINE.parents[1] / "nextflow.config"


def test_runtime_policy_is_fail_fast_and_reported():
    root = ROOT_CONFIG.read_text()
    pipeline = (PIPELINE / "nextflow.config").read_text()
    assert "nextflowVersion = '!>=26.04.6'" in root
    assert "shell         = ['/bin/bash', '-euo', 'pipefail']" in root
    assert "enabled = true" in root
    assert "task.exitStatus in [137, 140, 143]" in root
    assert "task.exitStatus in [137, 140, 143]" in pipeline
    assert "? 'retry' : 'terminate'" in root
    assert "? 'retry' : 'terminate'" in pipeline
    assert "errorStrategy 'ignore'" not in pipeline
    assert "withName: 'TAMA_COLLAPSE'" in pipeline
    assert "maxRetries = 3" in pipeline
    assert "128.GB * task.attempt" in pipeline
    assert "task.exitStatus in [137, 140, 143]" in pipeline


def test_all_pipeline_helpers_are_executable():
    helpers = (PIPELINE / "bin").glob("*.py")
    assert all(os.access(path, os.X_OK) for path in helpers)


def test_workflow_exposes_required_boundaries():
    main = (PIPELINE / "main.nf").read_text()
    subworkflows = "\n".join(path.read_text() for path in (PIPELINE / "subworkflows").glob("*.nf"))
    for name in (
        "PREPARE_LONG_READS", "ALIGN_LONG_READS", "COLLAPSE_LONG_READ_MODELS",
        "MERGE_LONG_READ_MODELS", "VALIDATE_COMBINED_MODELS", "RUN_DIAMOND_QC",
    ):
        assert f"{name}" in main or f"{name}" in subworkflows
    for name in (
        "VALIDATE_FASTQ", "RUN_PBCCS", "BAM_TO_FASTQ", "MINIMAP2_ALIGN",
        "SAMTOOLS_SORT_INDEX", "ALIGNMENT_QC", "VALIDATE_TAMA_OUTPUT",
    ):
        assert f"process {name}" in subworkflows or f"process {name}" in "\n".join(
            path.read_text() for path in (PIPELINE / "modules").glob("*.nf")
        )


def test_entrypoint_uses_schema_and_keeps_optional_outputs_guarded():
    main = (PIPELINE / "main.nf").read_text()
    align = (PIPELINE / "subworkflows" / "align_long_reads.nf").read_text()
    schema = (PIPELINE / "nextflow_schema.json").read_text()

    assert "validateParameters()" in main
    assert "INVENTORY_LONG_READS" in main
    assert '"$schema": "https://json-schema.org/draft/2020-12/schema"' in schema
    assert "build_versions = channel.empty()" in align
    assert "BUILD_MINIMAP2_INDEX.out.versions" not in align.split("emit:", 1)[1]
    assert "COLLECT_LONG_READ_SOFTWARE_VERSIONS" in main
