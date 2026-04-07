#!/usr/bin/env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
"""
Set up a Nextflow full-annotation experiment using the Ensembl assembly registry.

This is an adapter around start_pipeline_from_registry.py — it reuses all the
registry lookup, clade assignment, stable ID space allocation, and production-name
derivation logic from ensembl-genes, but instead of generating eHive configs it
produces:

  1. A Nextflow params JSON  (--params-file ready for our FullAnnotation_conf)
  2. An init_pipeline.pl command  (for the eHive wrapper FullAnnotation_conf.pm)
  3. A comparison manifest  (paths to reference GFF3 for parity testing)

Typical usage (on HPC with ENSCODE set):

  # Experiment on a single well-annotated bird genome
  python setup_experiment.py \\
      --gca GCA_003957565.2 \\
      --settings settings.json \\
      --outdir /hps/scratch/flicek/ensembl/genebuild/experiments/GCA_003957565.2 \\
      --nf-base-dir /nfs/production/flicek/ensembl/genebuild/ensembl-genes-nf/pipelines \\
      --nf-work-root /hps/scratch/flicek/ensembl/genebuild/nf_work \\
      --target-db-host mysql-ens-genebuild-prod \\
      --target-db-port 4527 \\
      --target-db-user ensrw \\
      --target-db-password SECRET \\
      --release 113

  # Dry-run (no registry connection — uses --init-file for metadata)
  python setup_experiment.py \\
      --gca GCA_003957565.2 \\
      --init-file taeniopygia_guttata.ini \\
      --outdir /scratch/experiment \\
      --nf-base-dir /path/to/ensembl-genes-nf/pipelines \\
      --nf-work-root /scratch/nf_work \\
      --dry-run

Requirements:
  - ENSCODE environment variable pointing to the parent of ensembl-genes/
  - Python packages: pymysql (for registry mode)
  - ensembl-genes must be on PYTHONPATH or installed

Environment variables used (same as start_pipeline_from_registry.py):
  GBS1/GBP1  — registry DB server (gb_assembly_metadata)
  GBS3/GBP3  — core DB server (target for gene loading)
  GBS4/GBP4  — pipeline DB server (eHive tracking)
  ENSCODE     — root of ensembl codebase checkouts
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Ensure ensembl-genes is importable (add to path if ENSCODE is set)
# ---------------------------------------------------------------------------
def _add_ensembl_genes_to_path() -> None:
    enscode = os.environ.get("ENSCODE")
    if enscode:
        genes_src = Path(enscode) / "ensembl-genes" / "src" / "python"
        if genes_src.exists():
            sys.path.insert(0, str(genes_src))


_add_ensembl_genes_to_path()

try:
    from ensembl.genes.info_from_registry.start_pipeline_from_registry import (
        add_generated_data,
        load_settings,
        get_server_settings_main,
        custom_loading,
    )
    from ensembl.genes.info_from_registry.assign_species_prefix import get_species_prefix
    from ensembl.genes.info_from_registry.assign_stable_space import get_stable_space
    from ensembl.genes.info_from_registry.taxonomy_helper import (
        assign_clade,
        assign_clade_info_custom_loading,
    )
    REGISTRY_AVAILABLE = True
except ImportError as e:
    REGISTRY_AVAILABLE = False
    logger.warning(
        f"ensembl-genes not importable ({e}). "
        "Registry lookup disabled — use --init-file or --metadata-json instead."
    )


# ---------------------------------------------------------------------------
# Clade → Nextflow param mapping
# ---------------------------------------------------------------------------
# These map clade_settings.json keys to the param names our FullAnnotation_conf expects.
# HPC-local paths come from clade_settings; we just pass them through.

VERTEBRATE_CLADES = {
    "primates", "rodentia", "mammalia", "marsupials",
    "aves", "lepidosauria", "testudines", "crocodylia",
    "amphibians", "teleostei", "sharks", "distant_vertebrate (chordata)",
}


def _uniprot_fasta_from_clade(clade: str, clade_settings: Dict, enscode: str) -> Optional[str]:
    """
    Resolve the UniProt protein FASTA path for a given clade.
    The main pipeline uses per-clade UniProt subsets stored on HPS.
    """
    # clade_settings already has 'protein_file' for non-vertebrate clades.
    # For vertebrate clades, the uniprot set is determined by uniprot_set key.
    protein_file = clade_settings.get("protein_file")
    if protein_file and Path(protein_file).exists():
        return protein_file

    uniprot_set = clade_settings.get("uniprot_set")
    if uniprot_set:
        # Vertebrate protein sets are in the genebuild_virtual_user area
        base = "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets"
        candidate = Path(base) / f"{uniprot_set}_uniprot_vertebrata_sp.fa"
        if candidate.exists():
            return str(candidate)

    return None


def _igtr_fasta_from_clade(clade_settings: Dict) -> Optional[str]:
    """Resolve IG/TR protein FASTA from clade settings."""
    ig_tr_file = clade_settings.get("ig_tr_fasta_file")
    if not ig_tr_file:
        return None
    base = "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets"
    full = Path(base) / ig_tr_file
    return str(full) if full.exists() else None


# ---------------------------------------------------------------------------
# Core adapter function: all_output_params → Nextflow params
# ---------------------------------------------------------------------------

def build_nextflow_params(
    info_dict: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """
    Translate the enriched registry metadata dict (from add_generated_data)
    into a Nextflow params dict for FullAnnotation_conf / main.nf.

    This is the key mapping layer — it replaces edit_config_main() and
    build_annotation_commands() from the original script.
    """
    gca = info_dict["assembly_accession"]
    assembly_name = info_dict.get("assembly_name", "").replace(" ", "_")
    species_name = info_dict.get("species_name", "")
    production_name = info_dict.get("production_name", species_name)

    outdir = Path(args.outdir)

    # ----------------------------------------------------------------
    # Core DB name — derived the same way as start_pipeline_from_registry.py
    # ----------------------------------------------------------------
    dbname_accession = gca.replace(".", "v").replace("_", "").lower()
    release = str(args.release)
    dbowner = args.dbowner
    target_db_name = f"{dbowner}_{dbname_accession}_core_{release}_1"

    # ----------------------------------------------------------------
    # Clade-specific params
    # ----------------------------------------------------------------
    clade = info_dict.get("clade", "")
    repbase_library = info_dict.get("repbase_library", "vertebrates")
    igtr_proteins = _igtr_fasta_from_clade(info_dict)
    enscode = os.environ.get("ENSCODE", "")
    uniprot_fasta = _uniprot_fasta_from_clade(clade, info_dict, enscode)

    # ----------------------------------------------------------------
    # RNA-seq and long-read sample sheets
    # (paths resolved by start_pipeline_from_registry via get_info_for_pipeline_main)
    # ----------------------------------------------------------------
    sample_sheet = info_dict.get("rnaseq_summary_file", "")
    long_read_sample_sheet = info_dict.get("long_read_summary_file", "")
    # Genus-level RNA-seq as fallback
    rnaseq_genus = info_dict.get("rnaseq_summary_file_genus", "")

    # ----------------------------------------------------------------
    # Projection source DB (from clade_settings)
    # ----------------------------------------------------------------
    projection_source = info_dict.get("projection_source_production_name", "homo_sapiens")

    # ----------------------------------------------------------------
    # Assemble the params dict
    # ----------------------------------------------------------------
    params = {
        # Assembly identity
        "assembly_accession":          gca,
        "assembly_name":               assembly_name,
        "assembly_refseq_accession":   info_dict.get("assembly_refseq_accession", ""),

        # Infrastructure
        "nf_base_dir":                 args.nf_base_dir,
        "nextflow_work_root":          args.nf_work_root,
        "outdir":                      str(outdir),

        # Repeat masking
        "repbase_library":             repbase_library,
        "custom_repeat_library":       info_dict.get("repeatmodeler_library", ""),
        "repeat_species":              clade or "vertebrates",

        # Evidence inputs
        "sample_sheet":                sample_sheet,
        "long_read_sample_sheet":      long_read_sample_sheet,
        "uniprot_fasta":               uniprot_fasta or "",
        "igtr_proteins":               igtr_proteins or "",

        # Projection
        "projection_source_production_name": projection_source,

        # Species metadata
        "species_name":                production_name,
        "stable_id_prefix":            info_dict.get("stable_id_prefix", ""),
        "stable_id_start":             info_dict.get("stable_id_start", 1),

        # Target core DB
        "target_db_host":              args.target_db_host or "",
        "target_db_port":              args.target_db_port,
        "target_db_user":              args.target_db_user or "",
        "target_db_password":          args.target_db_password or "",
        "target_db_name":              target_db_name,
        "assembly_version":            assembly_name,

        # Metadata
        "release":                     release,
        "clade":                       clade,
        "taxon_id":                    info_dict.get("taxon_id", ""),
        "busco_group":                 info_dict.get("busco_group", ""),
    }

    # Strip empty-string optional fields so Nextflow uses defaults
    params = {k: v for k, v in params.items() if v != ""}

    return params


def build_init_pipeline_command(
    params: Dict[str, Any],
    args: argparse.Namespace,
    enscode: str,
) -> str:
    """
    Generate the init_pipeline.pl command that would launch the eHive
    FullAnnotation_conf.pm wrapper (for comparison / hybrid mode).
    """
    conf = "Bio::EnsEMBL::Analysis::Hive::Config::Nextflow::FullAnnotation_conf"
    pipeline_db = (
        f"-host {args.pipeline_db_host or 'PIPELINE_DB_HOST'} "
        f"-port {args.pipeline_db_port or 4527} "
        f"-user {args.target_db_user or 'ensrw'} "
        f"-pass SECRET "
        f"-dbname {params.get('species_name', 'experiment')}_fullannot_pipe"
    )
    lines = [
        f"init_pipeline.pl {conf} \\",
        f'    -pipeline_db                "{pipeline_db}" \\',
    ]
    skip = {"release", "clade", "taxon_id", "busco_group", "nf_base_dir",
            "nextflow_work_root", "target_db_host", "target_db_port",
            "target_db_user", "target_db_password", "target_db_name",
            "assembly_version", "stable_id_start"}
    for k, v in params.items():
        if k in skip:
            continue
        lines.append(f"    -{k:<35} {v} \\")
    # Add infra params
    lines += [
        f"    -nf_base_dir                {params.get('nf_base_dir', '')} \\",
        f"    -nextflow_work_root         {params.get('nextflow_work_root', '')} \\",
        f"    -target_db_host             {args.target_db_host or ''} \\",
        f"    -target_db_port             {args.target_db_port} \\",
        f"    -target_db_user             {args.target_db_user or ''} \\",
        f"    -target_db_name             {params.get('target_db_name', '')} \\",
        f"    -assembly_version           {params.get('assembly_version', '')} \\",
        f"    -stable_id_start            {params.get('stable_id_start', 1)}",
    ]
    return "\n".join(lines)


def build_comparison_manifest(
    info_dict: Dict[str, Any],
    params: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """
    Write a comparison manifest pointing to:
    - Where the new pipeline will write its final GFF3
    - Where the existing Ensembl reference GFF3 can be found
      (FTP path, for download if needed)

    This feeds directly into tests/parity/metrics/compare.py.
    """
    gca = info_dict["assembly_accession"]
    species_name = info_dict.get("species_name", "")
    production_name = info_dict.get("production_name", species_name)
    assembly_name = info_dict.get("assembly_name", "").replace(" ", "_")

    outdir = Path(args.outdir)
    new_pipeline_gff3_dir = outdir / "finalise_geneset" / "finalise_geneset"

    release = str(args.release)
    ftp_base = "https://ftp.ensembl.org/pub"
    ftp_gff3 = (
        f"{ftp_base}/release-{release}/gff3/{production_name}/"
        f"{production_name.capitalize()}.{assembly_name}.{release}.gff3.gz"
    )

    return {
        "gca": gca,
        "production_name": production_name,
        "assembly_name": assembly_name,
        "release": release,
        "new_pipeline": {
            "gff3_dir": str(new_pipeline_gff3_dir),
            "description": "Output of our Nextflow FullAnnotation pipeline",
        },
        "reference": {
            "ftp_gff3": ftp_gff3,
            "description": "Published Ensembl GFF3 (reference for parity testing)",
            "download_cmd": (
                f"wget -O reference.gff3.gz '{ftp_gff3}' && "
                f"gunzip reference.gff3.gz"
            ),
        },
        "parity_test_cmd": (
            f"python3 tests/parity/run_parity_comparison.py "
            f"--ref reference.gff3 "
            f"--test {new_pipeline_gff3_dir}/*.gff3 "
            f"--output parity_results_{gca}.json"
        ),
    }


# ---------------------------------------------------------------------------
# Metadata loading (registry or offline)
# ---------------------------------------------------------------------------

def load_metadata_from_registry(
    gca: str,
    settings: Dict,
    server_info: Dict,
) -> Dict[str, Any]:
    """Load and enrich metadata from the Ensembl assembly registry."""
    if not REGISTRY_AVAILABLE:
        raise RuntimeError(
            "ensembl-genes not importable. Use --metadata-json for offline mode."
        )
    info_dict = add_generated_data(server_info, gca, settings)
    # Resolve stable ID space
    info_dict["stable_id_prefix"] = get_species_prefix(
        info_dict["taxon_id"], server_info
    )
    info_dict["stable_id_start"] = get_stable_space(
        info_dict["taxon_id"],
        gca,
        info_dict["assembly_id"],
        server_info,
    )
    return info_dict


def load_metadata_from_json(json_path: str) -> Dict[str, Any]:
    """Load pre-built metadata from a JSON file (offline mode)."""
    with open(json_path) as f:
        return json.load(f)


def load_metadata_from_init_file(
    init_file: str,
    gca: str,
) -> Dict[str, Any]:
    """
    Load metadata from a simple key=value INI file (same format as start_pipeline_from_registry.py).
    Requires ENSCODE to be set for clade assignment.
    """
    if not REGISTRY_AVAILABLE:
        raise RuntimeError("ensembl-genes not importable.")
    settings = {"init_file": init_file}
    info_dict = custom_loading(settings)
    info_dict = assign_clade_info_custom_loading(info_dict)
    info_dict.setdefault("assembly_accession", gca)
    info_dict.setdefault("stable_id_prefix", "")
    info_dict.setdefault("stable_id_start", 1)
    return info_dict


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Set up a Nextflow annotation experiment from the Ensembl registry.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Metadata source (one of these three)
    meta_group = parser.add_mutually_exclusive_group()
    meta_group.add_argument(
        "--settings-file", metavar="JSON",
        help="Pipeline settings JSON (same as start_pipeline_from_registry.py --settings_file).",
    )
    meta_group.add_argument(
        "--metadata-json", metavar="JSON",
        help="Pre-built metadata JSON (output of a previous run or manual assembly).",
    )
    meta_group.add_argument(
        "--init-file", metavar="INI",
        help="key=value INI file for offline/custom metadata (no registry connection).",
    )

    parser.add_argument("--gca", required=True, help="GCA accession e.g. GCA_003957565.2")

    # Output / infrastructure
    parser.add_argument("--outdir", required=True, help="Experiment output directory (NFS/HPS path)")
    parser.add_argument("--nf-base-dir", required=True,
                        help="Parent of ensembl-genes-nf/pipelines/ on the filesystem")
    parser.add_argument("--nf-work-root", required=True,
                        help="Nextflow work directory root (scratch space)")
    parser.add_argument("--release", type=int, default=113, help="Ensembl release number")
    parser.add_argument("--dbowner", default="ensembl", help="DB owner prefix for core DB name")

    # Target core DB
    parser.add_argument("--target-db-host", help="MySQL host for target core DB")
    parser.add_argument("--target-db-port", type=int, default=3306)
    parser.add_argument("--target-db-user", help="MySQL user for target core DB")
    parser.add_argument("--target-db-password", default="", help="MySQL password")

    # eHive pipeline DB (for init_pipeline.pl command only)
    parser.add_argument("--pipeline-db-host", help="eHive pipeline tracking DB host")
    parser.add_argument("--pipeline-db-port", type=int, default=4527)

    # Output files
    parser.add_argument("--params-out", default="nextflow_params.json",
                        help="Output path for Nextflow --params-file JSON")
    parser.add_argument("--comparison-out", default="comparison_manifest.json",
                        help="Output path for parity comparison manifest")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print params without writing files or contacting registry")

    args = parser.parse_args()

    # ----------------------------------------------------------------
    # Load metadata
    # ----------------------------------------------------------------
    if args.metadata_json:
        logger.info("Loading metadata from JSON: %s", args.metadata_json)
        info_dict = load_metadata_from_json(args.metadata_json)
        info_dict.setdefault("assembly_accession", args.gca)

    elif args.init_file:
        logger.info("Loading metadata from INI file: %s", args.init_file)
        info_dict = load_metadata_from_init_file(args.init_file, args.gca)

    elif args.settings_file:
        logger.info("Loading metadata from registry using settings: %s", args.settings_file)
        settings = load_settings(args.settings_file)
        server_info = {
            "registry": {
                "db_host": os.environ.get("GBS1"),
                "db_user": settings.get("user_r"),
                "db_port": str(os.environ.get("GBP1")),
                "db_name": "gb_assembly_metadata",
                "password": settings.get("password", ""),
            }
        }
        server_settings = get_server_settings_main(settings)
        server_info.update(server_settings)
        info_dict = load_metadata_from_registry(args.gca, settings, server_info)

    else:
        parser.error("One of --settings-file, --metadata-json, or --init-file is required.")
        return  # unreachable but keeps type checkers happy

    # ----------------------------------------------------------------
    # Build Nextflow params
    # ----------------------------------------------------------------
    params = build_nextflow_params(info_dict, args)

    # ----------------------------------------------------------------
    # Build supporting outputs
    # ----------------------------------------------------------------
    init_cmd = build_init_pipeline_command(
        params, args, enscode=os.environ.get("ENSCODE", "")
    )
    comparison = build_comparison_manifest(info_dict, params, args)

    # ----------------------------------------------------------------
    # Print / write
    # ----------------------------------------------------------------
    if args.dry_run:
        print("\n=== Nextflow params (--params-file) ===")
        print(json.dumps(params, indent=2, default=str))
        print("\n=== init_pipeline.pl command (eHive wrapper) ===")
        print(init_cmd)
        print("\n=== Comparison manifest ===")
        print(json.dumps(comparison, indent=2, default=str))
        return

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    params_path = outdir / args.params_out
    with open(params_path, "w") as f:
        json.dump(params, f, indent=2, default=str)
    logger.info("Wrote Nextflow params to: %s", params_path)

    comparison_path = outdir / args.comparison_out
    with open(comparison_path, "w") as f:
        json.dump(comparison, f, indent=2, default=str)
    logger.info("Wrote comparison manifest to: %s", comparison_path)

    init_cmd_path = outdir / "init_pipeline_cmd.sh"
    init_cmd_path.write_text(f"#!/bin/bash\n{init_cmd}\n")
    init_cmd_path.chmod(0o755)
    logger.info("Wrote init_pipeline.pl command to: %s", init_cmd_path)

    print(f"\n=== Experiment set up for {args.gca} ===")
    print(f"Nextflow params:     {params_path}")
    print(f"Comparison manifest: {comparison_path}")
    print(f"eHive init command:  {init_cmd_path}")
    print(f"\nTo run the Nextflow pipeline directly:")
    print(
        f"  nextflow run {args.nf_base_dir}/full_annotation/main.nf \\\n"
        f"    -profile cluster \\\n"
        f"    -params-file {params_path} \\\n"
        f"    -work-dir {args.nf_work_root}/{args.gca}/work"
    )
    print(f"\nTo download the reference GFF3 for parity testing:")
    print(f"  {comparison['reference']['download_cmd']}")
    print(f"\nTo run parity comparison (after pipeline completes):")
    print(f"  {comparison['parity_test_cmd']}")


if __name__ == "__main__":
    main()
