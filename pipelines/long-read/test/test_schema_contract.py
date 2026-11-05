import json
from pathlib import Path


ROOT = Path(__file__).parents[1]


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
