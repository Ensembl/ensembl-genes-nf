#!/usr/bin/env python3
"""
Load a final-geneset GFF3 into an Ensembl core MySQL database.

Creates gene/transcript/exon/translation rows in the core schema.
Designed to run after finalise_geneset produces the canonical gene set.

Schema compatibility: Ensembl schema >=109 (current production schema).

Usage:
    load_gff3_to_core.py \\
        --gff3          final.gff3 \\
        --host          mysql-host \\
        --port          3306 \\
        --user          ensadmin \\
        --password      secret \\
        --dbname        homo_sapiens_core_109_38 \\
        --analysis      genebuild \\
        --coord-system  chromosome \\
        --assembly      GRCh38 \\
        --species-id    1
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pymysql
import pymysql.cursors

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GFF3 data model
# ---------------------------------------------------------------------------

@dataclass
class GffExon:
    seq_region_name: str
    start: int          # 1-based, inclusive
    end: int
    strand: int         # 1 or -1

@dataclass
class GffCds:
    seq_region_name: str
    start: int
    end: int
    strand: int
    phase: int          # GFF3 phase (0,1,2): bases to skip at the START of this CDS feature

@dataclass
class GffTranscript:
    transcript_id: str
    gene_id: str
    seq_region_name: str
    start: int
    end: int
    strand: int
    biotype: str
    source: str = "ensembl"
    exons: List[GffExon] = field(default_factory=list)
    cds_features: List[GffCds] = field(default_factory=list)
    attributes: Dict[str, str] = field(default_factory=dict)

    @property
    def is_coding(self) -> bool:
        return bool(self.cds_features)

@dataclass
class GffGene:
    gene_id: str
    seq_region_name: str
    start: int
    end: int
    strand: int
    biotype: str
    name: str = ""
    source: str = "ensembl"
    transcripts: List[GffTranscript] = field(default_factory=list)
    attributes: Dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# GFF3 parser
# ---------------------------------------------------------------------------

def _parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, _, v = part.partition("=")
        attrs[k.strip()] = v.strip()
    return attrs


def parse_gff3(path: str) -> List[GffGene]:
    """Parse a GFF3 file into a list of GffGene objects."""
    genes: Dict[str, GffGene] = {}
    transcripts: Dict[str, GffTranscript] = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqname, source, feature, start, end, score, strand_char, phase_char, attr_str = cols
            start = int(start)
            end = int(end)
            strand = 1 if strand_char == "+" else -1
            attrs = _parse_attrs(attr_str)
            feature_id = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            biotype = attrs.get("biotype", "protein_coding")

            if feature == "gene":
                genes[feature_id] = GffGene(
                    gene_id=feature_id,
                    seq_region_name=seqname,
                    start=start,
                    end=end,
                    strand=strand,
                    biotype=biotype,
                    name=attrs.get("Name", feature_id),
                    source=source,
                    attributes=attrs,
                )

            elif feature in ("mRNA", "transcript", "lnc_RNA", "pseudogenic_transcript",
                             "processed_pseudogene", "pseudogene"):
                tx = GffTranscript(
                    transcript_id=feature_id,
                    gene_id=parent,
                    seq_region_name=seqname,
                    start=start,
                    end=end,
                    strand=strand,
                    biotype=biotype,
                    source=source,
                    attributes=attrs,
                )
                transcripts[feature_id] = tx

            elif feature == "exon":
                if parent in transcripts:
                    transcripts[parent].exons.append(GffExon(
                        seq_region_name=seqname,
                        start=start,
                        end=end,
                        strand=strand,
                    ))

            elif feature == "CDS":
                gff3_phase = int(phase_char) if phase_char not in (".", "") else 0
                if parent in transcripts:
                    transcripts[parent].cds_features.append(GffCds(
                        seq_region_name=seqname,
                        start=start,
                        end=end,
                        strand=strand,
                        phase=gff3_phase,
                    ))

    # Attach transcripts to genes
    for tx in transcripts.values():
        if tx.gene_id in genes:
            genes[tx.gene_id].transcripts.append(tx)
        else:
            # Orphan transcript — create synthetic gene
            g = GffGene(
                gene_id=tx.gene_id or tx.transcript_id,
                seq_region_name=tx.seq_region_name,
                start=tx.start,
                end=tx.end,
                strand=tx.strand,
                biotype=tx.biotype,
                name=tx.gene_id or tx.transcript_id,
                source=tx.source,
            )
            g.transcripts.append(tx)
            genes[g.gene_id] = g

    # Sort exons and CDS features by genomic position
    for tx in transcripts.values():
        tx.exons.sort(key=lambda e: e.start)
        tx.cds_features.sort(key=lambda c: c.start)

    return list(genes.values())


# ---------------------------------------------------------------------------
# Coordinate math for Ensembl schema
# ---------------------------------------------------------------------------

def compute_translation_coords(tx: GffTranscript) -> Optional[Tuple[
    GffExon, int, GffExon, int  # start_exon, seq_start, end_exon, seq_end
]]:
    """
    Compute Ensembl translation coordinates:
      start_exon: exon containing the CDS start codon
      seq_start:  1-based offset from the exon start (on + strand) or
                  exon end (on - strand) to the first CDS base
      end_exon:   exon containing the CDS stop codon (last base before stop)
      seq_end:    analogous offset for the last CDS base

    Returns None if no CDS features.
    """
    if not tx.cds_features or not tx.exons:
        return None

    cds = sorted(tx.cds_features, key=lambda c: c.start)
    cds_start = cds[0].start   # genomic start of first CDS feature
    cds_end = cds[-1].end      # genomic end of last CDS feature

    exons = sorted(tx.exons, key=lambda e: e.start)

    # Find start_exon: first exon that contains cds_start
    start_exon = next((e for e in exons if e.start <= cds_start <= e.end), None)
    # Find end_exon: first exon that contains cds_end
    end_exon = next((e for e in exons if e.start <= cds_end <= e.end), None)

    if not start_exon or not end_exon:
        log.warning(f"Could not find start/end exon for {tx.transcript_id}")
        return None

    if tx.strand == 1:
        seq_start = cds_start - start_exon.start + 1
        seq_end   = cds_end   - end_exon.start   + 1
    else:
        seq_start = end_exon.end   - cds_end   + 1
        seq_end   = start_exon.end - cds_start + 1
        # On minus strand, start_exon is actually the exon with the highest genomic coord
        start_exon, end_exon = end_exon, start_exon

    return start_exon, seq_start, end_exon, seq_end


def exon_phase_from_cds(tx: GffTranscript) -> Dict[Tuple[int,int], Tuple[int,int]]:
    """
    Return a mapping of (exon.start, exon.end) → (phase, end_phase) for each exon.
    phase / end_phase = -1 for fully non-coding exons.
    """
    result: Dict[Tuple[int,int], Tuple[int,int]] = {}
    exons = sorted(tx.exons, key=lambda e: e.start)
    for e in exons:
        result[(e.start, e.end)] = (-1, -1)

    if not tx.cds_features:
        return result

    cds = sorted(tx.cds_features, key=lambda c: c.start)
    if tx.strand == -1:
        cds = list(reversed(cds))

    # Walk CDS features in transcript order, tracking running CDS length mod 3
    cumulative = 0
    for c in cds:
        # GFF3 phase = bases to SKIP at the start of this feature on the + strand.
        # For Ensembl: exon phase = (3 - gff3_phase) % 3 at the CDS-proximal end.
        # The simpler derivation: phase = cumulative % 3.
        exon_start_key = (c.start, c.end)
        phase = cumulative % 3
        cds_len = c.end - c.start + 1
        cumulative += cds_len
        end_phase = cumulative % 3

        # Find which exon(s) this CDS overlaps
        for e in exons:
            if e.start <= c.start and e.end >= c.end:
                result[(e.start, e.end)] = (phase, end_phase)
                break

    return result


# ---------------------------------------------------------------------------
# Ensembl core DB loader
# ---------------------------------------------------------------------------

class EnsemblCoreLoader:
    def __init__(self, conn: pymysql.Connection, analysis_id: int, species_id: int,
                 coord_system_id: int, seq_region_map: Dict[str, int]):
        self.conn = conn
        self.analysis_id = analysis_id
        self.species_id = species_id
        self.coord_system_id = coord_system_id
        self.seq_region_map = seq_region_map  # name → seq_region_id
        self._exon_cache: Dict[Tuple, int] = {}  # (sr_id, start, end, strand, phase, end_phase) → exon_id

    def _sr_id(self, name: str) -> int:
        if name not in self.seq_region_map:
            raise ValueError(f"Seq region '{name}' not in seq_region_map. "
                             f"Known: {list(self.seq_region_map)[:5]}")
        return self.seq_region_map[name]

    def load_gene(self, gene: GffGene) -> int:
        sr_id = self._sr_id(gene.seq_region_name)
        with self.conn.cursor() as cur:
            cur.execute(
                """INSERT INTO gene
                   (seq_region_id, seq_region_start, seq_region_end,
                    seq_region_strand, analysis_id, source, biotype,
                    is_current, canonical_transcript_id)
                   VALUES (%s,%s,%s,%s,%s,%s,%s,1,0)""",
                (sr_id, gene.start, gene.end, gene.strand,
                 self.analysis_id, gene.source, gene.biotype),
            )
            gene_db_id = cur.lastrowid

        # Load transcripts and collect their DB IDs
        canonical_tx_id = None
        for tx in gene.transcripts:
            tx_db_id = self.load_transcript(tx, gene_db_id)
            # Canonical transcript: the one with is_canonical=1 attribute, else first
            if tx.attributes.get("canonical_transcript") == "1" or canonical_tx_id is None:
                canonical_tx_id = tx_db_id

        # Update canonical_transcript_id
        if canonical_tx_id:
            with self.conn.cursor() as cur:
                cur.execute(
                    "UPDATE gene SET canonical_transcript_id=%s WHERE gene_id=%s",
                    (canonical_tx_id, gene_db_id),
                )

        return gene_db_id

    def load_transcript(self, tx: GffTranscript, gene_db_id: int) -> int:
        sr_id = self._sr_id(tx.seq_region_name)
        with self.conn.cursor() as cur:
            cur.execute(
                """INSERT INTO transcript
                   (gene_id, seq_region_id, seq_region_start, seq_region_end,
                    seq_region_strand, analysis_id, source, biotype, is_current,
                    canonical_translation_id)
                   VALUES (%s,%s,%s,%s,%s,%s,%s,%s,1,NULL)""",
                (gene_db_id, sr_id, tx.start, tx.end, tx.strand,
                 self.analysis_id, tx.source, tx.biotype),
            )
            tx_db_id = cur.lastrowid

        # Load exons
        phase_map = exon_phase_from_cds(tx)
        exon_db_ids = []
        exons_sorted = sorted(tx.exons, key=lambda e: e.start,
                              reverse=(tx.strand == -1))

        for rank, exon in enumerate(exons_sorted, start=1):
            phase, end_phase = phase_map.get((exon.start, exon.end), (-1, -1))
            exon_db_id = self._get_or_create_exon(
                sr_id, exon.start, exon.end, tx.strand, phase, end_phase
            )
            exon_db_ids.append((rank, exon_db_id, exon))
            with self.conn.cursor() as cur:
                cur.execute(
                    """INSERT INTO exon_transcript (exon_id, transcript_id, rank)
                       VALUES (%s, %s, %s)""",
                    (exon_db_id, tx_db_id, rank),
                )

        # Load translation if coding
        if tx.is_coding:
            self._load_translation(tx, tx_db_id, exon_db_ids)

        return tx_db_id

    def _get_or_create_exon(self, sr_id, start, end, strand, phase, end_phase) -> int:
        key = (sr_id, start, end, strand, phase, end_phase)
        if key in self._exon_cache:
            return self._exon_cache[key]
        with self.conn.cursor() as cur:
            cur.execute(
                """INSERT INTO exon
                   (seq_region_id, seq_region_start, seq_region_end,
                    seq_region_strand, phase, end_phase, is_current, is_constitutive)
                   VALUES (%s,%s,%s,%s,%s,%s,1,0)""",
                (sr_id, start, end, strand, phase, end_phase),
            )
            exon_id = cur.lastrowid
        self._exon_cache[key] = exon_id
        return exon_id

    def _load_translation(self, tx: GffTranscript, tx_db_id: int,
                          exon_db_ids: List[Tuple[int, int, GffExon]]) -> Optional[int]:
        coords = compute_translation_coords(tx)
        if not coords:
            return None
        start_exon, seq_start, end_exon, seq_end = coords

        # Find DB IDs for start and end exons by coordinates
        start_exon_db_id = end_exon_db_id = None
        for _rank, exon_db_id, exon in exon_db_ids:
            if exon.start == start_exon.start and exon.end == start_exon.end:
                start_exon_db_id = exon_db_id
            if exon.start == end_exon.start and exon.end == end_exon.end:
                end_exon_db_id = exon_db_id

        if not start_exon_db_id or not end_exon_db_id:
            log.warning(f"Could not resolve translation exon IDs for {tx.transcript_id}")
            return None

        with self.conn.cursor() as cur:
            cur.execute(
                """INSERT INTO translation
                   (transcript_id, seq_start, start_exon_id, seq_end, end_exon_id)
                   VALUES (%s,%s,%s,%s,%s)""",
                (tx_db_id, seq_start, start_exon_db_id, seq_end, end_exon_db_id),
            )
            tl_id = cur.lastrowid

        # Set canonical_translation_id on transcript
        with self.conn.cursor() as cur:
            cur.execute(
                "UPDATE transcript SET canonical_translation_id=%s WHERE transcript_id=%s",
                (tl_id, tx_db_id),
            )

        return tl_id


# ---------------------------------------------------------------------------
# Database setup helpers
# ---------------------------------------------------------------------------

def ensure_coord_system(conn, species_id: int, cs_name: str, assembly: str) -> int:
    """Get or create a coord_system row; return coord_system_id."""
    with conn.cursor() as cur:
        cur.execute(
            "SELECT coord_system_id FROM coord_system WHERE name=%s AND version=%s AND species_id=%s",
            (cs_name, assembly, species_id),
        )
        row = cur.fetchone()
        if row:
            return row["coord_system_id"]
        cur.execute(
            """INSERT INTO coord_system (species_id, name, version, rank, attrib)
               VALUES (%s,%s,%s,1,'default_version,sequence_level')""",
            (species_id, cs_name, assembly),
        )
        return cur.lastrowid


def ensure_seq_regions(conn, coord_system_id: int, genome_fai: Optional[str] = None,
                       synonyms_tsv: Optional[str] = None) -> Dict[str, int]:
    """
    Load seq_region rows from a .fai index file and/or synonyms TSV.
    Returns mapping of seq_region name → seq_region_id.
    """
    sr_map: Dict[str, int] = {}

    # First, load any existing seq_regions for this coord_system
    with conn.cursor() as cur:
        cur.execute(
            "SELECT seq_region_id, name FROM seq_region WHERE coord_system_id=%s",
            (coord_system_id,),
        )
        for row in cur.fetchall():
            sr_map[row["name"]] = row["seq_region_id"]

    if genome_fai:
        with open(genome_fai) as fh:
            for line in fh:
                name, length, *_ = line.strip().split("\t")
                if name not in sr_map:
                    with conn.cursor() as cur:
                        cur.execute(
                            "INSERT INTO seq_region (name, coord_system_id, length) VALUES (%s,%s,%s)",
                            (name, coord_system_id, int(length)),
                        )
                        sr_map[name] = cur.lastrowid
                        log.debug(f"  Inserted seq_region {name} ({length} bp)")

    # Add synonym aliases so GFF3 names like 'chr1' resolve to 'Chr1' etc.
    if synonyms_tsv:
        with open(synonyms_tsv) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                canonical, synonym = parts[0], parts[1]
                if synonym not in sr_map and canonical in sr_map:
                    sr_map[synonym] = sr_map[canonical]

    return sr_map


def ensure_analysis(conn, logic_name: str) -> int:
    """Get or create an analysis row; return analysis_id."""
    with conn.cursor() as cur:
        cur.execute("SELECT analysis_id FROM analysis WHERE logic_name=%s", (logic_name,))
        row = cur.fetchone()
        if row:
            return row["analysis_id"]
        cur.execute(
            "INSERT INTO analysis (created, logic_name, module) VALUES (NOW(),%s,'Nextflow')",
            (logic_name,),
        )
        return cur.lastrowid


def ensure_meta(conn, species_id: int, assembly: str, species_name: str):
    """Write essential meta key/value pairs."""
    entries = [
        ("assembly.default", assembly),
        ("species.production_name", species_name),
        ("schema.version", "109"),
        ("genebuild.method", "full_genebuild_nextflow"),
        ("genebuild.start_date", __import__("datetime").date.today().isoformat()),
    ]
    with conn.cursor() as cur:
        for key, val in entries:
            cur.execute(
                "SELECT meta_id FROM meta WHERE meta_key=%s AND species_id=%s",
                (key, species_id),
            )
            if not cur.fetchone():
                cur.execute(
                    "INSERT INTO meta (species_id, meta_key, meta_value) VALUES (%s,%s,%s)",
                    (species_id, key, val),
                )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Load final-geneset GFF3 into Ensembl core DB.")
    ap.add_argument("--gff3",          required=True)
    ap.add_argument("--host",          required=True)
    ap.add_argument("--port",          type=int, default=3306)
    ap.add_argument("--user",          required=True)
    ap.add_argument("--password",      default="")
    ap.add_argument("--dbname",        required=True)
    ap.add_argument("--analysis",      default="ensembl_genebuild",
                    help="Logic name for the genebuild analysis entry")
    ap.add_argument("--coord-system",  default="chromosome")
    ap.add_argument("--assembly",      required=True,
                    help="Assembly version string, e.g. GRCh38")
    ap.add_argument("--species-id",    type=int, default=1)
    ap.add_argument("--species-name",  default="species",
                    help="Ensembl production name, e.g. homo_sapiens")
    ap.add_argument("--genome-fai",    default=None,
                    help="Path to samtools .fai index for seq_region loading")
    ap.add_argument("--synonyms-tsv",  default=None,
                    help="Seq-region synonyms TSV (from load_assembly)")
    ap.add_argument("--stats-json",    default="load_stats.json",
                    help="Output: JSON file with load statistics")
    return ap.parse_args()


def main():
    args = parse_args()
    log.info(f"Connecting to {args.host}:{args.port}/{args.dbname}")

    conn = pymysql.connect(
        host=args.host,
        port=args.port,
        user=args.user,
        password=args.password,
        database=args.dbname,
        cursorclass=pymysql.cursors.DictCursor,
        autocommit=False,
    )

    try:
        log.info("Setting up analysis, coord_system, seq_regions...")
        analysis_id = ensure_analysis(conn, args.analysis)
        cs_id = ensure_coord_system(conn, args.species_id, args.coord_system, args.assembly)
        sr_map = ensure_seq_regions(conn, cs_id, args.genome_fai, args.synonyms_tsv)
        ensure_meta(conn, args.species_id, args.assembly, args.species_name)

        log.info(f"Parsing GFF3: {args.gff3}")
        genes = parse_gff3(args.gff3)
        log.info(f"  {len(genes)} genes to load")

        loader = EnsemblCoreLoader(
            conn=conn,
            analysis_id=analysis_id,
            species_id=args.species_id,
            coord_system_id=cs_id,
            seq_region_map=sr_map,
        )

        stats = {"genes": 0, "transcripts": 0, "coding_transcripts": 0,
                 "exons": 0, "translations": 0, "skipped_genes": 0}

        for i, gene in enumerate(genes):
            if i % 1000 == 0 and i > 0:
                log.info(f"  Loaded {i}/{len(genes)} genes...")
                conn.commit()

            try:
                loader.load_gene(gene)
                stats["genes"] += 1
                for tx in gene.transcripts:
                    stats["transcripts"] += 1
                    if tx.is_coding:
                        stats["coding_transcripts"] += 1
            except Exception as e:
                log.warning(f"Skipping gene {gene.gene_id}: {e}")
                conn.rollback()
                stats["skipped_genes"] += 1
                continue

        conn.commit()

        # Count exons and translations loaded
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) AS n FROM exon")
            stats["exons"] = cur.fetchone()["n"]
            cur.execute("SELECT COUNT(*) AS n FROM translation")
            stats["translations"] = cur.fetchone()["n"]

        log.info(f"Load complete: {stats}")
        with open(args.stats_json, "w") as fh:
            json.dump(stats, fh, indent=2)

    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
