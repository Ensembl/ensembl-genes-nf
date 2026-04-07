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
Set up a parity experiment using the Ensembl assembly registry.

Architecture reminder
---------------------
Each STAGE of the annotation is a Nextflow pipeline (pipelines/load_assembly/,
pipelines/rnaseq/, pipelines/consolidate/, etc.).  The ORCHESTRATION of those
stages — running them in order, passing outputs between them, managing retries —
is the job of the eHive PipeConfig (FullAnnotation_conf.pm) now and will move
to Prefect later.

There is no monolithic "full_annotation/main.nf".  This script does NOT produce
a Nextflow --params-file.  It produces:

  1. init_pipeline_cmd.sh   — the init_pipeline.pl command for FullAnnotation_conf.pm
                              (runs the eHive DAG, which launches each NF sub-pipeline)
  2. prefect_params.json    — equivalent params dict in the shape the future Prefect
                              flow will expect (stub for now, tracks the migration target)
  3. comparison_manifest.json — FTP URL + expected output paths for parity testing
                              against the published Ensembl reference GFF3

Typical usage (on HPC with ENSCODE set):

  # With live registry connection
  python setup_experiment.py \\
      --gca GCA_003957565.2 \\
      --settings-file settings.json \\
      --outdir /hps/scratch/flicek/ensembl/genebuild/experiments/GCA_003957565.2 \\
      --nf-base-dir /nfs/production/flicek/ensembl/genebuild/ensembl-genes-nf/pipelines \\
      --nf-work-root /hps/scratch/flicek/ensembl/genebuild/nf_work \\
      --pipeline-db-host mysql-ens-genebuild-prod --pipeline-db-port 4527 \\
      --target-db-host mysql-ens-genebuild-prod   --target-db-port 4527 \\
      --target-db-user ensrw --release 113

  # Dry-run with pre-built metadata JSON (no registry connection)
  python setup_experiment.py \\
      --gca GCA_003957565.2 \\
      --metadata-json taeniopygia_guttata_meta.json \\
      --outdir /scratch/experiment \\
      --nf-base-dir /path/to/ensembl-genes-nf/pipelines \\
      --nf-work-root /scratch/nf_work \\
      --dry-run

Then run the experiment:

  # 1. Create an empty core DB schema
  mysql -u ensrw -p -e "CREATE DATABASE <target_db_name>;"
  mysql -u ensrw -p <target_db_name> < ensembl-core-schema.sql

  # 2. Init the eHive pipeline
  bash init_pipeline_cmd.sh

  # 3. Start workers
  beekeeper.pl -url $EHIVE_URL -loop

  # 4. After completion — download reference and compare
  bash comparison_manifest.json | jq -r '.reference.download_cmd' | bash
  python scripts/run_parity_comparison.py \\
      --ref reference.gff3 \\
      --test <outdir>/finalise_geneset/finalise_geneset/*.gff3 \\
      --output parity_GCA_003957565.2.json
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
# Ensure ensembl-genes is importable
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
        assign_clade_info_custom_loading,
    )
    REGISTRY_AVAILABLE = True
except ImportError as e:
    REGISTRY_AVAILABLE = False
    logger.warning(
        f"ensembl-genes not importable ({e}). "
        "Registry lookup disabled — use --metadata-json or --init-file."
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _igtr_fasta_from_clade(info_dict: Dict) -> Optional[str]:
    ig_tr_file = info_dict.get("ig_tr_fasta_file")
    if not ig_tr_file:
        return None
    base = "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets"
    return str(Path(base) / ig_tr_file)


def _protein_fasta_from_clade(info_dict: Dict) -> Optional[str]:
    """Return clade-specific UniProt FASTA, same logic as start_pipeline_from_registry."""
    # For non-vertebrate clades, it's in protein_file directly
    pf = info_dict.get("protein_file")
    if pf:
        return pf
    # For vertebrate clades — uniprot_set key selects the file
    uniprot_set = info_dict.get("uniprot_set")
    if uniprot_set:
        base = "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets"
        return str(Path(base) / f"{uniprot_set}_uniprot_vertebrata_sp.fa")
    return None


# ---------------------------------------------------------------------------
# eHive init_pipeline.pl command
# ---------------------------------------------------------------------------

def build_init_pipeline_command(
    info_dict: Dict[str, Any],
    args: argparse.Namespace,
) -> str:
    """
    Generate the init_pipeline.pl command for FullAnnotation_conf.pm.

    This is the primary output of this script.  Running it initialises the
    eHive tracking database and loads the DAG.  Workers launched with
    beekeeper.pl will then:
      1. Run each Nextflow sub-pipeline (load_assembly, repeat_masking, etc.)
      2. Pass GFF3 paths between them via channel-2 dataflow
      3. Load the final geneset into the target core DB

    The per-pipeline Nextflow params (genome_fasta, sample_sheet, etc.) are
    resolved inside FullAnnotation_conf.pm from the assembly_accession +
    assembly_name, using the same publishDir conventions as each NF pipeline.
    """
    gca = info_dict["assembly_accession"]
    assembly_name = info_dict.get("assembly_name", "").replace(" ", "_")
    species_name = info_dict.get("species_name", "")
    production_name = info_dict.get("production_name", species_name)
    clade = info_dict.get("clade", "vertebrates")

    dbname_accession = gca.replace(".", "v").replace("_", "").lower()
    release = str(args.release)
    dbowner = args.dbowner
    target_db_name = f"{dbowner}_{dbname_accession}_core_{release}_1"
    pipe_db_name = f"{dbowner}_{dbname_accession}_fullannot_pipe_{release}"

    enscode = os.environ.get("ENSCODE", "ENSCODE_NOT_SET")
    conf = "Bio::EnsEMBL::Analysis::Hive::Config::Nextflow::FullAnnotation_conf"

    pipeline_db = (
        f"-host {args.pipeline_db_host} "
        f"-port {args.pipeline_db_port} "
        f"-user {args.target_db_user} "
        f"-pass {args.target_db_password or 'SECRET'} "
        f"-dbname {pipe_db_name}"
    )

    # Optional evidence — only include params that are resolvable
    optional_params = {}
    protein_fasta = _protein_fasta_from_clade(info_dict)
    if protein_fasta:
        optional_params["uniprot_fasta"] = protein_fasta
    igtr = _igtr_fasta_from_clade(info_dict)
    if igtr:
        optional_params["igtr_proteins"] = igtr
    if info_dict.get("rnaseq_summary_file"):
        optional_params["sample_sheet"] = info_dict["rnaseq_summary_file"]
    if info_dict.get("long_read_summary_file"):
        optional_params["long_read_sample_sheet"] = info_dict["long_read_summary_file"]
    if info_dict.get("repeatmodeler_library"):
        optional_params["custom_repeat_library"] = info_dict["repeatmodeler_library"]
    if info_dict.get("projection_source_production_name"):
        optional_params["projection_source_production_name"] = \
            info_dict["projection_source_production_name"]
    if info_dict.get("assembly_refseq_accession"):
        optional_params["assembly_refseq_accession"] = info_dict["assembly_refseq_accession"]

    lines = [
        f"init_pipeline.pl {conf} \\",
        f'    -pipeline_db         "{pipeline_db}" \\',
        f"    -assembly_accession  {gca} \\",
        f"    -assembly_name       {assembly_name} \\",
        f"    -repbase_library     {info_dict.get('repbase_library', 'vertebrates')} \\",
        f"    -repeat_species      {clade} \\",
        f"    -nf_base_dir         {args.nf_base_dir} \\",
        f"    -nextflow_work_root  {args.nf_work_root} \\",
        f"    -outdir              {args.outdir} \\",
        f"    -target_db_host      {args.target_db_host} \\",
        f"    -target_db_port      {args.target_db_port} \\",
        f"    -target_db_user      {args.target_db_user} \\",
        f"    -target_db_name      {target_db_name} \\",
        f"    -assembly_version    {assembly_name} \\",
        f"    -species_name        {production_name} \\",
        f"    -stable_id_prefix    '{info_dict.get('stable_id_prefix', '')}' \\",
        f"    -stable_id_start     {info_dict.get('stable_id_start', 1)} \\",
        f"    -release             {release} \\",
    ]
    for k, v in optional_params.items():
        lines.append(f"    -{k:<28} {v} \\")

    # Remove trailing backslash from last line
    lines[-1] = lines[-1].rstrip(" \\")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Prefect params (future orchestration target)
# ---------------------------------------------------------------------------

def build_prefect_params(
    info_dict: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """
    Produce the equivalent parameter dict for the future Prefect flow.

    The Prefect flow will replace eHive as the orchestrator: instead of
    init_pipeline.pl + beekeeper.pl, a Prefect flow will directly call
    nextflow.run() for each sub-pipeline in the right order.

    The parameter structure mirrors FullAnnotation_conf.pm's default_options
    so the migration is a near-direct translation.  Stub for now — the shape
    is what matters for planning the Prefect implementation.
    """
    gca = info_dict["assembly_accession"]
    assembly_name = info_dict.get("assembly_name", "").replace(" ", "_")
    production_name = info_dict.get("production_name", info_dict.get("species_name", ""))
    dbname_accession = gca.replace(".", "v").replace("_", "").lower()
    release = str(args.release)
    target_db_name = f"{args.dbowner}_{dbname_accession}_core_{release}_1"

    return {
        # Identifies the experiment
        "gca": gca,
        "assembly_name": assembly_name,
        "assembly_refseq_accession": info_dict.get("assembly_refseq_accession", ""),
        "species_name": production_name,
        "clade": info_dict.get("clade", ""),
        "taxon_id": info_dict.get("taxon_id", ""),
        "busco_group": info_dict.get("busco_group", ""),
        "stable_id_prefix": info_dict.get("stable_id_prefix", ""),
        "stable_id_start": info_dict.get("stable_id_start", 1),
        "release": release,

        # Infrastructure — Prefect flow will pass these into each NF sub-pipeline
        "nf_base_dir": args.nf_base_dir,
        "nf_work_root": args.nf_work_root,
        "outdir": args.outdir,

        # Per-stage NF params that differ by species (resolved from registry)
        "repbase_library": info_dict.get("repbase_library", "vertebrates"),
        "custom_repeat_library": info_dict.get("repeatmodeler_library", ""),
        "uniprot_fasta": _protein_fasta_from_clade(info_dict) or "",
        "igtr_proteins": _igtr_fasta_from_clade(info_dict) or "",
        "sample_sheet": info_dict.get("rnaseq_summary_file", ""),
        "long_read_sample_sheet": info_dict.get("long_read_summary_file", ""),
        "projection_source_production_name": info_dict.get(
            "projection_source_production_name", "homo_sapiens"
        ),

        # Target core DB (for gff3_to_core stage)
        "target_db_host": args.target_db_host,
        "target_db_port": args.target_db_port,
        "target_db_user": args.target_db_user,
        "target_db_name": target_db_name,
        "assembly_version": assembly_name,

        # NOTE: when Prefect replaces eHive, the flow will:
        #   1. Call load_assembly NF pipeline  → resolves genome_fasta + synonyms_tsv
        #   2. Call repeat_masking NF pipeline → resolves softmasked_fasta
        #   3. Fan out to all annotation NF sub-pipelines in parallel
        #   4. Call consolidate NF pipeline
        #   5. Call utr_addition NF pipeline
        #   6. Call finalise_geneset NF pipeline
        #   7. Call gff3_to_core NF pipeline
        # The paths between stages are passed as Prefect task outputs,
        # replacing eHive's channel-2 dataflow mechanism.
        "_migration_note": (
            "Prefect tasks will replace HiveRunNextflow calls. "
            "Each task returns the output GFF3 path and passes it to the next task. "
            "Replace beekeeper.pl with prefect run full_annotation_flow.py"
        ),
    }


# ---------------------------------------------------------------------------
# Comparison manifest
# ---------------------------------------------------------------------------

def build_comparison_manifest(
    info_dict: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """
    Write a manifest for running the parity comparison after the experiment.

    Points to:
    - Where our pipeline will write the final GFF3 (based on publishDir conventions)
    - The FTP URL for the published reference Ensembl GFF3
    - The exact run_parity_comparison.py command to use
    """
    gca = info_dict["assembly_accession"]
    production_name = info_dict.get("production_name", "")
    assembly_name = info_dict.get("assembly_name", "").replace(" ", "_")
    release = str(args.release)

    # Our pipeline publishes the final GFF3 here (from finalise_geneset publishDir)
    outdir = Path(args.outdir)
    our_gff3_dir = outdir / "finalise_geneset" / "finalise_geneset"

    ftp_base = "https://ftp.ensembl.org/pub"
    species_cap = production_name.split("_")[0].capitalize() + "_" + "_".join(production_name.split("_")[1:]) if production_name else ""
    ftp_gff3 = (
        f"{ftp_base}/release-{release}/gff3/{production_name}/"
        f"{species_cap}.{assembly_name}.{release}.gff3.gz"
    )

    return {
        "gca": gca,
        "production_name": production_name,
        "assembly_name": assembly_name,
        "release": release,
        "our_pipeline_gff3_dir": str(our_gff3_dir),
        "reference": {
            "ftp_gff3": ftp_gff3,
            "description": "Published Ensembl GFF3 — download this to use as parity reference",
            "download_cmd": (
                f"wget -O reference_{gca}.gff3.gz '{ftp_gff3}' && "
                f"gunzip reference_{gca}.gff3.gz"
            ),
        },
        "parity_cmd": (
            f"python scripts/run_parity_comparison.py \\\n"
            f"    --ref reference_{gca}.gff3 \\\n"
            f"    --test {our_gff3_dir}/*.gff3 \\\n"
            f"    --label 'Nextflow pipeline vs Ensembl release {release}' \\\n"
            f"    --output parity_results_{gca}.json"
        ),
    }


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def load_metadata_from_registry(
    gca: str, settings: Dict, server_info: Dict,
) -> Dict[str, Any]:
    if not REGISTRY_AVAILABLE:
        raise RuntimeError("ensembl-genes not importable. Use --metadata-json.")
    info_dict = add_generated_data(server_info, gca, settings)
    info_dict["stable_id_prefix"] = get_species_prefix(
        info_dict["taxon_id"], server_info
    )
    info_dict["stable_id_start"] = get_stable_space(
        info_dict["taxon_id"], gca, info_dict["assembly_id"], server_info,
    )
    return info_dict


def load_metadata_from_json(json_path: str, gca: str) -> Dict[str, Any]:
    with open(json_path) as f:
        d = json.load(f)
    d.setdefault("assembly_accession", gca)
    d.setdefault("stable_id_prefix", "")
    d.setdefault("stable_id_start", 1)
    return d


def load_metadata_from_init_file(init_file: str, gca: str) -> Dict[str, Any]:
    if not REGISTRY_AVAILABLE:
        raise RuntimeError("ensembl-genes not importable.")
    info_dict = custom_loading({"init_file": init_file})
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
        description=(
            "Set up an annotation parity experiment. "
            "Outputs the init_pipeline.pl command for FullAnnotation_conf.pm "
            "(eHive orchestrates the individual Nextflow sub-pipelines) "
            "and the equivalent Prefect params for the future migration."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Metadata source
    meta_group = parser.add_mutually_exclusive_group()
    meta_group.add_argument(
        "--settings-file", metavar="JSON",
        help="Pipeline settings JSON for live registry connection.",
    )
    meta_group.add_argument(
        "--metadata-json", metavar="JSON",
        help="Pre-built metadata JSON (from a previous run or manual assembly).",
    )
    meta_group.add_argument(
        "--init-file", metavar="INI",
        help="key=value INI for offline metadata (no registry connection).",
    )

    parser.add_argument("--gca", required=True, help="GCA accession e.g. GCA_003957565.2")

    # Infrastructure
    parser.add_argument("--outdir", required=True,
                        help="Root output dir for this experiment (NFS/HPS path)")
    parser.add_argument("--nf-base-dir", required=True,
                        help="Path to ensembl-genes-nf/pipelines/ on the filesystem")
    parser.add_argument("--nf-work-root", required=True,
                        help="Nextflow work directory root (scratch space)")
    parser.add_argument("--release", type=int, default=113)
    parser.add_argument("--dbowner", default="ensembl")

    # eHive pipeline DB
    parser.add_argument("--pipeline-db-host", default="",
                        help="MySQL host for the eHive tracking DB")
    parser.add_argument("--pipeline-db-port", type=int, default=4527)

    # Target core DB (for gff3_to_core loading)
    parser.add_argument("--target-db-host", default="",
                        help="MySQL host for the target Ensembl core DB")
    parser.add_argument("--target-db-port", type=int, default=3306)
    parser.add_argument("--target-db-user", default="ensrw")
    parser.add_argument("--target-db-password", default="")

    # Output filenames
    parser.add_argument("--prefect-params-out", default="prefect_params.json")
    parser.add_argument("--comparison-out", default="comparison_manifest.json")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print outputs without writing files or contacting registry")

    args = parser.parse_args()

    # ----------------------------------------------------------------
    # Load metadata
    # ----------------------------------------------------------------
    if args.metadata_json:
        info_dict = load_metadata_from_json(args.metadata_json, args.gca)
    elif args.init_file:
        info_dict = load_metadata_from_init_file(args.init_file, args.gca)
    elif args.settings_file:
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
        server_info.update(get_server_settings_main(settings))
        info_dict = load_metadata_from_registry(args.gca, settings, server_info)
    else:
        parser.error("One of --settings-file, --metadata-json, or --init-file is required.")
        return

    # ----------------------------------------------------------------
    # Build outputs
    # ----------------------------------------------------------------
    init_cmd = build_init_pipeline_command(info_dict, args)
    prefect_params = build_prefect_params(info_dict, args)
    comparison = build_comparison_manifest(info_dict, args)

    # ----------------------------------------------------------------
    # Print / write
    # ----------------------------------------------------------------
    if args.dry_run:
        print("\n=== eHive init_pipeline.pl command ===")
        print(init_cmd)
        print("\n=== Prefect params (future orchestration) ===")
        print(json.dumps(prefect_params, indent=2, default=str))
        print("\n=== Comparison manifest ===")
        print(json.dumps(comparison, indent=2, default=str))
        return

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    init_cmd_path = outdir / "init_pipeline_cmd.sh"
    init_cmd_path.write_text(f"#!/bin/bash\n{init_cmd}\n")
    init_cmd_path.chmod(0o755)
    logger.info("Wrote eHive init command to: %s", init_cmd_path)

    prefect_path = outdir / args.prefect_params_out
    with open(prefect_path, "w") as f:
        json.dump(prefect_params, f, indent=2, default=str)
    logger.info("Wrote Prefect params to: %s", prefect_path)

    comparison_path = outdir / args.comparison_out
    with open(comparison_path, "w") as f:
        json.dump(comparison, f, indent=2, default=str)
    logger.info("Wrote comparison manifest to: %s", comparison_path)

    gca = args.gca
    print(f"\n=== Experiment set up for {gca} ===")
    print(f"\nStep 1 — create an empty core DB schema:")
    print(f"  mysql -u {args.target_db_user} -p -e \\")
    print(f"    \"CREATE DATABASE {prefect_params['target_db_name']};\"")
    print(f"  mysql -u {args.target_db_user} -p {prefect_params['target_db_name']} \\")
    print(f"    < /path/to/ensembl-core-schema.sql")
    print(f"\nStep 2 — init the eHive DAG (this does NOT run any compute yet):")
    print(f"  bash {init_cmd_path}")
    print(f"\nStep 3 — start eHive workers (these launch the NF sub-pipelines):")
    print(f"  beekeeper.pl -url $EHIVE_URL -loop")
    print(f"\nStep 4 — after pipeline completes, download reference and compare:")
    print(f"  {comparison['reference']['download_cmd']}")
    print(f"  {comparison['parity_cmd']}")
    print(f"\nPrefect params (for future migration): {prefect_path}")


if __name__ == "__main__":
    main()
