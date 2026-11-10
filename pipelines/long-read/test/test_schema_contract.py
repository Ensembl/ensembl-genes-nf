import json
import re
from pathlib import Path


ROOT = Path(__file__).parents[1]


CONCRETE_VERSION_PROCESSES = {
    "ALIGN_FASTQ",
    "BAM_TO_FASTQ",
    "BUILD_MINIMAP_INDEX",
    "COUNT_CONTIG_SHARD_COVERAGE",
    "COUNT_SOURCE_COVERAGE",
    "DIAMOND_BLASTP",
    "FLAIR_TRANSCRIPTOME",
    "INDEX_CONTIG_SHARD",
    "INDEX_SORTED_ALIGNMENT",
    "ISOQUANT_ANNOTATION_FREE",
    "ISOQUANT_REFERENCE_GUIDED",
    "PREP_DIAMOND_DB",
    "RUN_PBCCS",
    "SAMTOOLS_FLAGSTAT",
    "SAMTOOLS_IDXSTATS",
    "SAMTOOLS_QC_SUMMARY",
    "SAMTOOLS_VIEW_ALIGNMENT",
    "SAMTOOLS_VIEW_SHARD",
    "SORT_FASTQ_ALIGNMENT",
    "SORT_INPUT_BAM",
    "SPLIT_BAM_BY_CONTIG",
    "STRINGTIE2_COLLAPSE",
    "STRINGTIE2_MERGE",
    "STRINGTIE3_COLLAPSE",
    "STRINGTIE3_MERGE",
    "TMERGE_COLLAPSE",
    "TMERGE_MERGE",
    "TRANSGRAM_ANNOTATION_FREE",
    "TRANSGRAM_ANNOTATION_GUIDED",
    "VALIDATE_BAM_INPUT",
    "VALIDATE_CONTIG_SHARDS",
    "VALIDATE_SORTED_ALIGNMENT",
}


def test_versions_only_report_concrete_tools():
    process_pattern = re.compile(r"(?ms)^process\s+([A-Z][A-Z0-9_]*)\s*\{.*?^}\s*$")
    version_output = re.compile(r"path ['\"]versions\.yml['\"], emit: versions")
    forbidden_label = re.compile(
        r"(?:pipeline|python|pinned-runtime|stub|explicit|prepared_inputs|sha256sum|longest_atg)"
    )

    for module in (ROOT / "modules").glob("*.nf"):
        source = module.read_text()
        for match in process_pattern.finditer(source):
            process_name, block = match.group(1), match.group(0)
            if not version_output.search(block):
                continue
            assert process_name in CONCRETE_VERSION_PROCESSES
            version_lines = [line for line in block.splitlines() if "versions.yml" in line]
            assert not any(forbidden_label.search(line) for line in version_lines)


def test_schemas_are_json_and_backend_set_is_unique():
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    assert schema["properties"]["model_backends"]["uniqueItems"] is True
    assert "string" in schema["properties"]["model_backends"]["type"]
    assert "string" in schema["properties"]["backend_scopes"]["type"]
    assert schema["properties"]["fastq_cache_dir"]["format"] == "directory-path"
    assert schema["properties"]["tama_max_retries"]["maximum"] == 5
    mode_rule = schema["allOf"][0]["if"]["properties"]["mode"]
    assert mode_rule.get("enum") == ["fastq"] or mode_rule.get("const") == "fastq"
    assert schema["properties"]["preparation_stages"]["items"]["enum"]
    assert "reference_metadata" in schema["properties"]
    assert set(schema["properties"]["tool_containers"]["required"]) == {
        "minimap2", "samtools", "tama", "tmerge", "pbccs", "stringtie2", "stringtie3",
        "isoquant", "flair", "bambu", "transgram", "diamond", "python_helper"
    }
    assert "transgram" in schema["properties"]["model_backends"]["items"]["enum"]
    assert schema["properties"]["transgram_data_type"]["enum"] == ["auto", "alpha", "beta"]
    assert "sha256" in schema["properties"]["tool_containers"]["properties"]["minimap2"]["pattern"]


def test_retrospective_requires_manifest():
    schema = json.loads((ROOT / "retrospective_schema.json").read_text())
    assert "retrospective_manifest" in schema["required"]
    assert schema["properties"]["model_backends"]["items"]["enum"]
    assert "transgram" in schema["properties"]["model_backends"]["items"]["enum"]
