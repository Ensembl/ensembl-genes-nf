from pathlib import Path
import os


PIPELINE = Path(__file__).parents[1]
ROOT_CONFIG = PIPELINE.parents[1] / "nextflow.config"


def test_runtime_policy_and_schema_are_present():
    root = ROOT_CONFIG.read_text()
    pipeline = (PIPELINE / "nextflow.config").read_text()
    assert "nextflowVersion = '!>=26.04.6'" in root
    assert "shell         = ['/bin/bash', '-euo', 'pipefail']" in root
    assert "enabled = true" in root
    assert "task.attempt <= 5 ? 'retry' : 'ignore'" in pipeline
    assert "shard_mode = 'contig'" in pipeline
    assert '"default": "contig"' in (PIPELINE / "nextflow_schema.json").read_text()
    assert "withName: '.*'" in pipeline
    assert "withName: 'TAMA_COLLAPSE'" in pipeline
    assert "maxRetries = 5" in pipeline
    assert "tama_memory_small = '128.GB'" in pipeline
    assert "errorStrategy = { task.attempt <= 5 ? 'retry' : 'ignore' }" in pipeline


def test_all_pipeline_helpers_are_executable():
    helpers = (PIPELINE / "bin").glob("*.py")
    assert all(os.access(path, os.X_OK) for path in helpers)
    assert os.access(PIPELINE / "bin" / "bam_to_alignment_gtf.sh", os.X_OK)


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

    alignment = (PIPELINE / "modules" / "minimap2_align.nf").read_text()
    align_workflow = (PIPELINE / "subworkflows" / "align_long_reads.nf").read_text()
    assert "docker://niemasd/minimap2_samtools:2.28_1.20" in alignment
    assert "| samtools sort" in alignment
    assert "alignment.sam" not in alignment
    assert "MINIMAP2_ALIGN.out.bam" in align_workflow
    assert "SAMTOOLS_SORT_INDEX(MINIMAP2_ALIGN.out.sam)" not in align_workflow


def test_contig_manifest_stays_keyed_to_its_shard_directory():
    split = (PIPELINE / "modules" / "split_contigs.nf").read_text()
    collapse = (PIPELINE / "subworkflows" / "collapse_long_read_models.nf").read_text()
    assert "path('contig_manifest.tsv'), emit: shards" in split
    assert ".out.shards.flatMap" in collapse
    assert ".out.shards.combine" not in collapse


def test_tama_shards_can_be_skipped_without_dropping_the_accession():
    collapse = (PIPELINE / "subworkflows" / "collapse_long_read_models.nf").read_text()
    schema = (PIPELINE / "nextflow_schema.json").read_text()
    assert '"skip_tama_shards"' in schema
    assert "skip_tama_shards" in collapse
    assert '"skip_model_shards"' in schema
    assert "backend_skip_keys('tama'" in collapse
    assert "TAMA_COLLAPSE(tama_bams, reference)" in collapse
    assert "shard_bams" in collapse


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
    assert "community.wave.seqera.io/library/pip_pyfaidx_setuptools_six_pruned:b4fcc6bbecfa2bea" in tmerge
    assert "bed12_to_gtf.py" in tmerge and "gtf_to_bed12.py" in tmerge
    assert "--input tmerge_input.gtf --output tmerge_output.gtf" in tmerge
    assert "--tmPrefix" not in tmerge


def test_model_backends_run_from_split_bams():
    config = (PIPELINE / "nextflow.config").read_text()
    schema = (PIPELINE / "nextflow_schema.json").read_text()
    collapse = (PIPELINE / "subworkflows" / "collapse_long_read_models.nf").read_text()
    assert "model_backend = 'tama'" in config
    assert '"model_backend"' in schema
    for backend in ("tama", "stringtie2", "stringtie3", "tmerge", "all"):
        assert backend in schema
    assert "STRINGTIE2_COLLAPSE(stringtie2_bams)" in collapse
    assert "STRINGTIE3_COLLAPSE(stringtie3_bams)" in collapse
    assert "BAM_TO_ALIGNMENT_GTF(tmerge_bams)" in collapse
    assert "TMERGE_COLLAPSE(tmerge_gtf)" in collapse
    assert "tuple('tama', meta, shard, bed)" in collapse
    assert "tuple('stringtie2', meta, shard, bed)" in collapse
    assert "tuple('stringtie3', meta, shard, bed)" in collapse
    assert "tuple('tmerge', meta, shard, bed)" in collapse
    assert "VALIDATE_BACKEND_BED(raw_beds)" in collapse
    tmerge = (PIPELINE / "modules" / "tmerge_collapse.nf").read_text()
    stringtie2 = (PIPELINE / "modules" / "stringtie2_collapse.nf").read_text()
    bam_to_gtf = (PIPELINE / "modules" / "bam_to_alignment_gtf.nf").read_text()
    assert "path(bam), path(bai)" in bam_to_gtf
    assert "process BAM_TO_ALIGNMENT_GTF" in bam_to_gtf
    assert "depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1" in bam_to_gtf
    assert "bam_to_alignment_gtf.sh" in bam_to_gtf
    converter = (PIPELINE / "bin" / "bam_to_alignment_gtf.sh").read_text()
    assert 'gene_id \\\"' in converter
    assert 'gene_id \\\"" read_id' in converter
    assert 'transcript_id \\\"" read_id' in converter
    assert converter.index('gene_id \\\"" read_id') < converter.index('transcript_id \\\"" read_id')
    assert "path(reads_gtf)" in tmerge
    assert "bam_to_alignment_gtf.py" not in tmerge
    assert "--input ${reads_gtf} --output ${prefix}.gtf" in tmerge
    assert "--tmPrefix" not in tmerge
    assert "depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0" in stringtie2
    assert "community.wave.seqera.io/library/stringtie" not in stringtie2


def test_all_backends_have_independent_merge_and_finalisation_contracts():
    main = (PIPELINE / "main.nf").read_text()
    merge = (PIPELINE / "subworkflows" / "merge_long_read_models.nf").read_text()
    validate = (PIPELINE / "subworkflows" / "validate_combined_models.nf").read_text()
    assert "MERGE_LONG_READ_MODELS(COLLAPSE_LONG_READ_MODELS.out.beds)" in main
    assert '"${backend}@@${meta.id}"' in merge
    assert '"${backend}@@${params.cohort_id}"' in merge
    assert "tuple val(backend), val(accession), path(beds" in (PIPELINE / "modules" / "tama_merge_accession.nf").read_text()
    assert "tuple val(backend), val(cohort_id), path(beds" in (PIPELINE / "modules" / "tama_merge.nf").read_text()
    assert "CANONICALISE_COMBINED_MODELS(canonical_input)" in validate
    assert "VALIDATE_LONG_READ_MODELS" in validate
    canonical = (PIPELINE / "modules" / "canonicalise_models.nf").read_text()
    assert '"${backend}_combined_models.bed"' in canonical
    assert '"${backend}_combined_models.sha256"' in canonical


def test_diamond_qc_receives_each_finalised_backend():
    main = (PIPELINE / "main.nf").read_text()
    diamond = (PIPELINE / "subworkflows" / "run_diamond_qc.nf").read_text()
    assert "if (run_diamond_validation)" in main
    assert "qc_bed = combined_bed.map { _backend, meta, bed -> tuple(meta, bed) }" in diamond
    assert "EXTRACT_COMBINED_TRANSCRIPTS(qc_bed, reference)" in diamond


def test_backend_bed_validation_writes_status_without_shell_awk():
    module = (PIPELINE / "modules" / "validate_backend_bed.nf").read_text()
    validator = (PIPELINE / "bin" / "validate_tama_bed.py").read_text()
    assert "--backend '${backend}' --status-output backend_status.tsv" in module
    assert "awk" not in module
    assert 'p.add_argument("--status-output")' in validator
    assert 'Path(a.status_output).write_text' in validator


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
