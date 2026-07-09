"""Persist and restore complete pipeline state across stage processes.

Each stage CLI runs in its own process. To stay byte-identical to the in-memory
monolith, a stage must see exactly the state the corresponding monolith method
would have seen on ``self``. Rather than reimplement any stage logic, every stage
loads the full :class:`~hpp.pipeline.PangenomeMappingPipeline` state from a
directory, runs the unmodified monolith method, and writes the updated state back
out. Nextflow passes the state directory from one process to the next.

The state directory holds the named artifacts from DEVELOPMENT.md §1
(``original_genes.jsonl``, ``syntenic_blocks.jsonl``, ``mapped_genes.jsonl`` …)
plus the carry-over reports and the resolved run config. Files are written only
when the corresponding attribute is populated, and loaded only if present, so a
stage transparently carries forward whatever earlier stages produced.
"""

from dataclasses import fields
from pathlib import Path
from typing import Optional, Set

from hpp.config import MappingConfig
from hpp.pipeline import PangenomeMappingPipeline
from hpp.serialize import (
    load_syntenic_map,
    read_json,
    read_jsonl,
    read_keyed_jsonl,
    save_syntenic_map,
    write_json,
    write_jsonl,
    write_keyed_jsonl,
    encode,
    decode,
)

CONFIG_FILE = "run_config.json"

# Fields of MappingConfig that need special (de)serialization.
_PATH_FIELDS = {"ref_fasta", "ref_gff", "target_fasta", "output_gff",
                "output_stats", "output_report", "temp_dir"}
_SET_FIELDS = {"chromosomes", "biotypes"}


# --- config ---------------------------------------------------------------

def config_to_dict(config: MappingConfig) -> dict:
    out = {}
    for f in fields(config):
        value = getattr(config, f.name)
        if value is None:
            out[f.name] = None
        elif f.name in _PATH_FIELDS:
            out[f.name] = str(value)
        elif f.name in _SET_FIELDS:
            out[f.name] = sorted(value)
        else:
            out[f.name] = value
    return out


def config_from_dict(data: dict) -> MappingConfig:
    kwargs = dict(data)
    for name in _SET_FIELDS:
        if kwargs.get(name) is not None:
            kwargs[name] = set(kwargs[name])
    # Paths are reconverted from str by MappingConfig.__post_init__.
    return MappingConfig(**kwargs)


def save_config(outdir: Path, config: MappingConfig) -> None:
    write_json(Path(outdir) / CONFIG_FILE, config_to_dict(config))


def load_config(indir: Path) -> MappingConfig:
    return config_from_dict(read_json(Path(indir) / CONFIG_FILE))


# --- full pipeline state --------------------------------------------------

# (attribute name, filename, kind) where kind is one of:
#   "keyed"   -> Dict[str, dataclass]  via keyed jsonl
#   "list"    -> List[dataclass]       via jsonl
#   "synmap"  -> SyntenicMap           via blocks jsonl (+ build_index on load)
#   "json"    -> plain dict / dataclass-free object
#   "obj"     -> single dataclass (encode/decode) e.g. SexChromosomeMap
_STATE_SPEC = [
    ("original_genes", "original_genes.jsonl", "keyed"),
    ("syntenic_map", "syntenic_blocks.jsonl", "synmap"),
    ("sex_chrom_map", "sex_chrom_map.json", "obj"),
    ("mapped_genes", "mapped_genes.jsonl", "keyed"),
    ("mapping_results", "mapping_results.jsonl", "list"),
    ("paralog_resolution_report", "paralog_report.json", "json"),
    ("validation_results", "validation_results.jsonl", "keyed"),
    ("initial_validation_results", "initial_validation_results.jsonl", "keyed"),
    ("protein_qc_results", "protein_qc.jsonl", "keyed"),
    ("protein_qc_report", "protein_qc_report.json", "json"),
    ("pre_refinement_genes", "pre_refinement_genes.jsonl", "keyed"),
    ("refinement_report", "refinement_report.json", "json"),
    ("audit_traces", "audit_traces.json", "json"),
]


def _is_empty(value) -> bool:
    return value is None or (isinstance(value, (dict, list)) and len(value) == 0)


def dump_state(pipeline: PangenomeMappingPipeline, outdir: Path) -> None:
    """Write every populated state attribute of ``pipeline`` to ``outdir``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_config(outdir, pipeline.config)
    for attr, fname, kind in _STATE_SPEC:
        value = getattr(pipeline, attr, None)
        if _is_empty(value):
            continue
        path = outdir / fname
        if kind == "keyed":
            write_keyed_jsonl(path, value)
        elif kind == "list":
            write_jsonl(path, value)
        elif kind == "synmap":
            save_syntenic_map(path, value)
        elif kind == "obj":
            write_json(path, encode(value))
        elif kind == "json":
            write_json(path, value)


def load_state(pipeline: PangenomeMappingPipeline, indir: Path) -> None:
    """Hydrate ``pipeline`` attributes from any state files present in ``indir``."""
    indir = Path(indir)
    for attr, fname, kind in _STATE_SPEC:
        path = indir / fname
        if not path.exists():
            continue
        if kind == "keyed":
            setattr(pipeline, attr, read_keyed_jsonl(path))
        elif kind == "list":
            setattr(pipeline, attr, read_jsonl(path))
        elif kind == "synmap":
            setattr(pipeline, attr, load_syntenic_map(path))
        elif kind == "obj":
            setattr(pipeline, attr, decode(read_json(path)))
        elif kind == "json":
            setattr(pipeline, attr, read_json(path))


def hydrate_pipeline(indir: Optional[Path]) -> PangenomeMappingPipeline:
    """Build a pipeline from a state dir's config and hydrate its attributes."""
    indir = Path(indir)
    config = load_config(indir)
    pipeline = PangenomeMappingPipeline(config)
    load_state(pipeline, indir)
    return pipeline
