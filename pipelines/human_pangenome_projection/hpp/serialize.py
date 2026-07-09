"""Lossless serialization of pipeline state to/from JSON-safe structures.

The staged pipeline persists every inter-stage object to disk so that each stage
can reconstruct its inputs without rerunning earlier stages. To produce output
byte-identical to the in-memory monolith, this serialization must round-trip the
relevant state exactly.

Design: a single generic encoder/decoder driven by a type registry. Every
dataclass is encoded as ``{"__type__": ClassName, <fields...>}`` and rebuilt via
``Cls(**fields)``. The ``Strand`` enum, ``set`` and ``tuple`` get explicit tags
because JSON has no native representation. Plain dicts (e.g. a feature's
``attributes``) are passed through as ordinary JSON objects.

Deliberately excluded fields (derived caches / not consumed downstream, so
excluding them keeps artifacts compact without changing any final output):
  * ``SyntenicBlock.cached_offset_map`` and ``SyntenicMap._*`` indexes — rebuilt
    by :meth:`SyntenicMap.build_index` on load.
  * ``FeatureMappingResult.source_blocks`` — never read after projection (only
    statistics counts consume mapping results); rebuilt as ``[]``.
"""

import json
from dataclasses import fields, is_dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from hpp.models import (
    CDS,
    Exon,
    Gene,
    GenomicInterval,
    Strand,
    SyntenicBlock,
    SyntenicMap,
    Transcript,
    UTR,
)
from hpp.stages.mapping import (
    FeatureMappingResult,
    GeneMappingResult,
    TranscriptMappingResult,
)
from hpp.stages.sex_chromosome import SexChromosomeMap
from hpp.validation.structural import (
    CodonValidation,
    GeneValidation,
    SpliceSiteValidation,
    TranscriptValidation,
)
from hpp.validation.protein import GeneProteinQC, TranscriptProteinQC

# Registry of every dataclass that may cross a stage boundary.
_REGISTRY = {
    cls.__name__: cls
    for cls in (
        GenomicInterval,
        Exon,
        CDS,
        UTR,
        Transcript,
        Gene,
        SyntenicBlock,
        SyntenicMap,
        FeatureMappingResult,
        TranscriptMappingResult,
        GeneMappingResult,
        SexChromosomeMap,
        SpliceSiteValidation,
        CodonValidation,
        TranscriptValidation,
        GeneValidation,
        TranscriptProteinQC,
        GeneProteinQC,
    )
}

# Per-class fields excluded from serialization (derived caches / unused).
_EXCLUDED_FIELDS = {
    "SyntenicBlock": {"cached_offset_map"},
    "FeatureMappingResult": {"source_blocks"},
}


def encode(obj: Any) -> Any:
    """Recursively convert a pipeline object into a JSON-safe structure."""
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, Strand):
        return {"__enum__": "Strand", "value": obj.value}
    if is_dataclass(obj) and not isinstance(obj, type):
        name = type(obj).__name__
        excluded = _EXCLUDED_FIELDS.get(name, set())
        out: Dict[str, Any] = {"__type__": name}
        for f in fields(obj):
            if f.name in excluded or f.name.startswith("_"):
                continue
            out[f.name] = encode(getattr(obj, f.name))
        return out
    if isinstance(obj, set):
        return {"__set__": [encode(x) for x in obj]}
    if isinstance(obj, tuple):
        return {"__tuple__": [encode(x) for x in obj]}
    if isinstance(obj, list):
        return [encode(x) for x in obj]
    if isinstance(obj, dict):
        # Plain dict (e.g. feature attributes): keys assumed to be strings.
        return {k: encode(v) for k, v in obj.items()}
    raise TypeError(f"Cannot encode object of type {type(obj)!r}")


def decode(obj: Any) -> Any:
    """Inverse of :func:`encode`."""
    if isinstance(obj, list):
        return [decode(x) for x in obj]
    if isinstance(obj, dict):
        if "__enum__" in obj:
            return Strand(obj["value"])
        if "__set__" in obj:
            return {decode(x) for x in obj["__set__"]}
        if "__tuple__" in obj:
            return tuple(decode(x) for x in obj["__tuple__"])
        if "__type__" in obj:
            cls = _REGISTRY[obj["__type__"]]
            kwargs = {k: decode(v) for k, v in obj.items() if k != "__type__"}
            return cls(**kwargs)
        return {k: decode(v) for k, v in obj.items()}
    return obj


# ---------------------------------------------------------------------------
# JSONL helpers for collections
# ---------------------------------------------------------------------------

def write_jsonl(path: Path, objects: Iterable[Any]) -> None:
    """Write an iterable of objects, one encoded JSON record per line."""
    path = Path(path)
    with path.open("w") as fh:
        for obj in objects:
            fh.write(json.dumps(encode(obj), separators=(",", ":")))
            fh.write("\n")


def read_jsonl(path: Path) -> List[Any]:
    """Read a JSONL file produced by :func:`write_jsonl`."""
    path = Path(path)
    out: List[Any] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(decode(json.loads(line)))
    return out


def write_keyed_jsonl(path: Path, mapping: Dict[str, Any]) -> None:
    """Write a ``Dict[str, obj]`` preserving key + insertion order."""
    path = Path(path)
    with path.open("w") as fh:
        for key, value in mapping.items():
            record = {"key": key, "value": encode(value)}
            fh.write(json.dumps(record, separators=(",", ":")))
            fh.write("\n")


def read_keyed_jsonl(path: Path) -> Dict[str, Any]:
    """Read a keyed JSONL file back into an insertion-ordered dict."""
    path = Path(path)
    out: Dict[str, Any] = {}
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            out[record["key"]] = decode(record["value"])
    return out


# ---------------------------------------------------------------------------
# Typed convenience loaders/savers
# ---------------------------------------------------------------------------

def save_syntenic_map(path: Path, smap: SyntenicMap) -> None:
    """Persist a SyntenicMap as JSONL of its blocks (indexes are derived)."""
    write_jsonl(path, smap.blocks)


def load_syntenic_map(path: Path) -> SyntenicMap:
    """Reconstruct a SyntenicMap and rebuild its spatial index + cs caches."""
    blocks = read_jsonl(path)
    smap = SyntenicMap(blocks=blocks)
    smap.build_index()
    return smap


def write_json(path: Path, obj: Any) -> None:
    """Write a plain report dict (already JSON-safe) with stable key order."""
    path = Path(path)
    with path.open("w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=False)


def read_json(path: Path) -> Any:
    path = Path(path)
    with Path(path).open() as fh:
        return json.load(fh)
