from pathlib import Path
import os


PIPELINE = Path(__file__).parents[1]
ROOT_CONFIG = PIPELINE.parents[1] / "nextflow.config"


def test_runtime_policy_is_bounded_and_ignores_exhausted_resource_kills():
    root = ROOT_CONFIG.read_text()
    pipeline = (PIPELINE / "nextflow.config").read_text()
    assert "nextflowVersion = '!>=26.04.6'" in root
    assert "shell         = ['/bin/bash', '-euo', 'pipefail']" in root
    assert "enabled = true" in root
    assert "task.exitStatus in [137, 140, 143]" in root
    assert "task.attempt <= 5 ? 'retry' : 'ignore'" in pipeline
    assert "? 'retry' : 'terminate'" in root
    assert "task.attempt <= 5 ? 'retry' : 'ignore'" in pipeline
    assert "shard_mode = 'contig'" in pipeline
    assert '"default": "contig"' in (PIPELINE / "nextflow_schema.json").read_text()
    assert "withName: '.*'" in pipeline
    assert "withName: 'TAMA_COLLAPSE'" in pipeline
    assert "maxRetries = 5" in pipeline
    assert "params.tama_memory_small" in pipeline
    assert "errorStrategy = { task.attempt <= 5 ? 'retry' : 'ignore' }" in pipeline


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


def test_contig_manifest_stays_keyed_to_its_shard_directory():
    split = (PIPELINE / "modules" / "split_contigs.nf").read_text()
    collapse = (PIPELINE / "subworkflows" / "collapse_long_read_models.nf").read_text()
    assert "path('contig_manifest.tsv'), emit: shards" in split
    assert ".out.shards.flatMap" in collapse
    assert ".out.shards.combine" not in collapse


def test_tama_soft_failures_have_explicit_status_and_diagnostics():
    tama = (PIPELINE / "modules" / "tama_collapse.nf").read_text()
    runner = (PIPELINE / "bin" / "run_tama_collapse.sh").read_text()
    assert "path('tama_status.tsv'), emit: status" in tama
    assert "path('tama_collapse.stderr'), emit: stderr" in tama
    assert "run_tama_collapse.sh" in tama
    assert "TAMA_FAILED" in runner
    assert "exit 0" in runner
    assert "|| return $?" in runner
    assert "IndexError: list index out of range" in runner
    assert "samtools view -bh -F 2308" in runner


def test_production_failures_are_non_terminal():
    config = (PIPELINE / "nextflow.config").read_text()
    assert "withName: '.*'" in config
    assert "withName: 'TAMA_COLLAPSE'" in config
    assert "errorStrategy = 'ignore'" in config
    assert "? 'retry' : 'terminate'" not in config


def test_tmerge_is_an_explicit_optional_merge_backend():
    config = (PIPELINE / "nextflow.config").read_text()
    schema = (PIPELINE / "nextflow_schema.json").read_text()
    tama = (PIPELINE / "modules" / "tama_merge.nf").read_text()
    tmerge = (PIPELINE / "modules" / "tmerge.nf").read_text()
    assert "merge_tool = 'tama'" in config
    assert '"merge_tool"' in schema and '"tama", "tmerge"' in schema
    assert "process TAMA_MERGE" in tama
    assert "params.merge_tool" not in tama
    assert "process TMERGE" in tmerge
    assert "community.wave.seqera.io/library/pip_tmerge:6cf60ff0bf166552" in tmerge
    assert "bed12_to_gtf.py" in tmerge and "gtf_to_bed12.py" in tmerge


def test_model_backends_run_from_split_bams():
    config = (PIPELINE / "nextflow.config").read_text()
    schema = (PIPELINE / "nextflow_schema.json").read_text()
    collapse = (PIPELINE / "subworkflows" / "collapse_long_read_models.nf").read_text()
    assert "model_backend = 'tama'" in config
    assert '"model_backend"' in schema
    for backend in ("tama", "stringtie2", "stringtie3", "tmerge", "all"):
        assert backend in schema
    assert "STRINGTIE2_COLLAPSE(contig_bams)" in collapse
    assert "STRINGTIE3_COLLAPSE(contig_bams)" in collapse
    assert "TMERGE_COLLAPSE(contig_bams)" in collapse
    tmerge = (PIPELINE / "modules" / "tmerge_collapse.nf").read_text()
    stringtie2 = (PIPELINE / "modules" / "stringtie2_collapse.nf").read_text()
    assert "path(bam), path(bai)" in tmerge
    assert "depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0" in stringtie2
    assert "community.wave.seqera.io/library/stringtie" not in stringtie2


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
