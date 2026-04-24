#!/usr/bin/env python3
"""
Set up a full-annotation experiment for the ensembl-genes-nf master pipeline.

Generates a ready-to-run shell script (`run_<assembly>.sh`) containing the
`nextflow run` command with all parameters resolved:
  - Queries ENA for short-read RNA-seq by taxon ID (or uses a provided sample sheet)
  - Derives layer priorities from enabled stages
  - Sets Augustus species model to closest trained model for the taxon
  - Optionally writes a Prefect migration JSON for future orchestration

Usage:
  python scripts/setup_experiment.py \\
    --accession      GCA_964261345.1 \\
    --assembly-name  mHetGlaV3 \\
    --outdir         /hps/scratch/flicek/ensembl/genebuild/jackt/hetgla \\
    --taxon-id       10181 \\
    --stable-id-prefix ENSHETG \\
    --stable-id-start  1
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Augustus species models — closest available trained model per taxon/clade
# ---------------------------------------------------------------------------
AUGUSTUS_SPECIES_MAP: dict[int, str] = {
    # Primates
    9606:  'human', 9544: 'human', 9598: 'human',
    # Rodents
    10090: 'mus_musculus',
    10116: 'rattus_norvegicus',
    10181: 'mus_musculus',   # Heterocephalus glaber
    9989:  'mus_musculus',   # Rodentia (clade)
    # Other mammals
    9913:  'cow',
    9796:  'horse',
    9823:  'pig',
    9940:  'sheep',
    9615:  'dog',
    40674: 'human',          # Mammalia (clade fallback)
    # Birds
    9031:  'chicken', 8782: 'chicken',
    # Fish
    7955:  'zebrafish', 7898: 'zebrafish',
    # Invertebrates
    7227:  'fly',
    6239:  'caenorhabditis',
    # Plants
    3702:  'arabidopsis', 4530: 'rice', 4577: 'maize',
    # Fungi
    5141:  'coprinus',
}
AUGUSTUS_DEFAULT = 'human'


def get_augustus_species(taxon_id: int | None) -> str:
    if taxon_id and taxon_id in AUGUSTUS_SPECIES_MAP:
        return AUGUSTUS_SPECIES_MAP[taxon_id]
    return AUGUSTUS_DEFAULT


# ---------------------------------------------------------------------------
# ENA RNA-seq lookup (best-effort)
# ---------------------------------------------------------------------------
def find_rnaseq_bioproject(taxon_id: int) -> str | None:
    import urllib.request
    query  = f"tax_eq({taxon_id})"
    fields = "study_accession,instrument_platform,library_strategy"
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search"
        f"?display=report&query={query}&domain=read&result=read_run"
        f"&fields={fields}&limit=10&format=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read())
        illumina = [r for r in data
                    if r.get("instrument_platform") == "ILLUMINA"
                    and r.get("library_strategy") in ("RNA-Seq", "EST")]
        if illumina:
            return illumina[0]["study_accession"]
    except Exception as exc:
        print(f"  [WARN] ENA lookup failed: {exc}", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# Layer priority map
# ---------------------------------------------------------------------------
def build_layer_priorities(stages: dict) -> dict:
    prio: dict[str, int] = {}
    if stages.get("long_read"):   prio["long_read"]          = 0
    if stages.get("targeted"):    prio["best_targeted"]       = 1
    if stages.get("rnaseq"):      prio["rnaseq"]              = 2
    if stages.get("projection"):  prio["projection"]          = 3
    if stages.get("refseq"):      prio["refseq_import"]       = 4
    prio["ab_initio"]             = 5
    if stages.get("genblast"):    prio["genblast_homology"]   = 6
    if stages.get("ncrna"):       prio["short_ncrna"]         = 7
    if stages.get("igtr"):        prio["igtr"]                = 7
    return prio


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Assembly
    parser.add_argument("--accession",     required=True, help="GCA accession")
    parser.add_argument("--assembly-name", required=True, help="Assembly name (no spaces)")
    parser.add_argument("--outdir",        required=True, help="Root output dir on HPC scratch")
    parser.add_argument("--nf-work-root",  default=None,  help="NF work dir (default: outdir/.nf_work)")
    parser.add_argument("--repo-dir",      default=".",   help="Path to ensembl-genes-nf repo")

    # Taxonomy
    parser.add_argument("--taxon-id",      type=int, default=None,
                        help="NCBI taxon ID — used for Augustus model + ENA RNA-seq lookup")
    parser.add_argument("--species-name",  default=None,  help="Scientific name (for core DB loading)")

    # RNA-seq (mutually exclusive)
    g_rna = parser.add_mutually_exclusive_group()
    g_rna.add_argument("--sample-sheet",          default=None, help="Local CSV sample sheet")
    g_rna.add_argument("--rnaseq-bioproject",      default=None, help="ENA BioProject accession")
    g_rna.add_argument("--rnaseq-run-accessions",  default=None, help="Comma-separated SRR/ERR accessions")

    # Protein homology (mutually exclusive)
    g_prot = parser.add_mutually_exclusive_group()
    g_prot.add_argument("--uniprot-fasta",     default=None, help="Local UniProt FASTA")
    g_prot.add_argument("--uniprot-taxon-id",  type=int, default=None, help="Taxon ID → auto-fetch from UniProt")

    # Other evidence
    parser.add_argument("--refseq-accession",   default=None, help="GCF accession for RefSeq import")
    parser.add_argument("--repbase-library",    default="vertebrates")
    parser.add_argument("--skip-repeatmodeler", action="store_true")
    parser.add_argument("--augustus-species",   default=None, help="Override Augustus species model")

    # Stable IDs / core DB
    parser.add_argument("--stable-id-prefix",  default="")
    parser.add_argument("--stable-id-start",   type=int, default=1)
    parser.add_argument("--db-host",            default=None)
    parser.add_argument("--db-user",            default=None)
    parser.add_argument("--db-name",            default=None)

    # Output
    parser.add_argument("--profile",        default="cluster", help="Nextflow profile (cluster/local/stub)")
    parser.add_argument("--output-script",  default=None, help="Output .sh path (default: run_<assembly>.sh)")

    args = parser.parse_args()

    outdir       = args.outdir.rstrip("/")
    nf_work_root = args.nf_work_root or f"{outdir}/.nf_work"
    repo_dir     = str(Path(args.repo_dir).resolve())
    script_path  = args.output_script or f"run_{args.assembly_name}.sh"

    # ── Resolve Augustus species ──────────────────────────────────────────
    augustus_species = args.augustus_species or get_augustus_species(args.taxon_id)
    print(f"Augustus species model: {augustus_species}", file=sys.stderr)

    # ── Resolve RNA-seq ────────────────────────────────────────────────────
    rnaseq_nf_arg = ""
    has_rnaseq    = False

    if args.sample_sheet:
        rnaseq_nf_arg = f"  --sample_sheet              {args.sample_sheet}"
        has_rnaseq    = True
    elif args.rnaseq_bioproject:
        rnaseq_nf_arg = f"  --rnaseq_bioproject         {args.rnaseq_bioproject}"
        has_rnaseq    = True
    elif args.rnaseq_run_accessions:
        rnaseq_nf_arg = f"  --rnaseq_run_accessions     {args.rnaseq_run_accessions}"
        has_rnaseq    = True
    elif args.taxon_id:
        print(f"No RNA-seq specified — querying ENA for taxon {args.taxon_id}...", file=sys.stderr)
        bioproject = find_rnaseq_bioproject(args.taxon_id)
        if bioproject:
            print(f"  Found: {bioproject}", file=sys.stderr)
            rnaseq_nf_arg = f"  --rnaseq_bioproject         {bioproject}"
            has_rnaseq    = True
        else:
            print("  [WARN] No RNA-seq found in ENA — running ab initio only.", file=sys.stderr)

    # ── Resolve UniProt ────────────────────────────────────────────────────
    uniprot_nf_arg = ""
    has_genblast   = False
    if args.uniprot_fasta:
        uniprot_nf_arg = f"  --uniprot_fasta             {args.uniprot_fasta}"
        has_genblast   = True
    elif args.uniprot_taxon_id:
        uniprot_nf_arg = f"  --uniprot_taxon_id          {args.uniprot_taxon_id}"
        has_genblast   = True

    # ── Layer priorities ───────────────────────────────────────────────────
    stages = {
        "rnaseq":   has_rnaseq,
        "genblast": has_genblast,
        "refseq":   bool(args.refseq_accession),
    }
    layer_prios_json = json.dumps(build_layer_priorities(stages))

    # ── Build nextflow run command ─────────────────────────────────────────
    nf_lines = [
        f"nextflow run {repo_dir}",
        f"  --assembly_accession        {args.accession}",
        f"  --assembly_name             {args.assembly_name}",
        f"  --outdir                    {outdir}",
        f"  --nf_work_root              {nf_work_root}",
        f"  --ab_initio_species         {augustus_species}",
        f"  --repbase_library           {args.repbase_library}",
        f"  --layer_priorities          '{layer_prios_json}'",
        f"  --stable_id_prefix          '{args.stable_id_prefix}'",
        f"  --stable_id_start           {args.stable_id_start}",
    ]

    if rnaseq_nf_arg:
        nf_lines.append(rnaseq_nf_arg)

    if uniprot_nf_arg:
        nf_lines.append(uniprot_nf_arg)

    if args.refseq_accession:
        nf_lines.append(f"  --assembly_refseq_accession {args.refseq_accession}")

    if args.skip_repeatmodeler:
        nf_lines.append("  --skip_repeatmodeler        true")

    if args.species_name:
        nf_lines.append(f"  --species_name              '{args.species_name}'")

    if args.db_host and args.db_name and args.db_user:
        nf_lines += [
            f"  --db_host                   {args.db_host}",
            f"  --db_user                   {args.db_user}",
            f"  --db_name                   {args.db_name}",
        ]

    nf_lines.append(f"  -profile                    {args.profile}")
    nf_lines.append("  -resume")

    nf_cmd = " \\\n".join(nf_lines)

    # ── Write run script ───────────────────────────────────────────────────
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    script_content = f"""#!/bin/bash
# ensembl-genes-nf full annotation
# Generated: {now}
# Assembly:  {args.accession} ({args.assembly_name})
# Taxon:     {args.taxon_id or 'not specified'}
#
# Launch inside tmux on the cluster:
#   tmux new -s {args.assembly_name.lower()}
#   salloc --ntasks=1 --mem=4G --time=7-00:00:00 --partition=standard
#   bash {script_path}

set -euo pipefail

echo "Starting annotation for {args.accession} at $(date)"
mkdir -p {outdir}

{nf_cmd}

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully at $(date)"
    echo "Final annotation: {outdir}/finalise_geneset/finalise_geneset/*.canonical.gff3"
else
    echo "Pipeline FAILED with exit code $EXIT_CODE at $(date)" >&2
    exit $EXIT_CODE
fi
"""

    Path(script_path).write_text(script_content)
    os.chmod(script_path, 0o755)

    # ── Prefect migration stub (future) ────────────────────────────────────
    prefect_params = {
        "_migration_note": (
            "This JSON will drive a Prefect flow once the eHive → Prefect migration "
            "is complete. Each key maps to a Prefect task input. The 'layer_priorities' "
            "dict maps to a task that calls the consolidate sub-pipeline."
        ),
        "assembly_accession":        args.accession,
        "assembly_name":             args.assembly_name,
        "outdir":                    outdir,
        "ab_initio_species":         augustus_species,
        "repbase_library":           args.repbase_library,
        "layer_priorities":          build_layer_priorities(stages),
        "rnaseq_source":             rnaseq_nf_arg.strip() if rnaseq_nf_arg else None,
        "uniprot_source":            uniprot_nf_arg.strip() if uniprot_nf_arg else None,
        "assembly_refseq_accession": args.refseq_accession,
        "stable_id_prefix":          args.stable_id_prefix,
        "stable_id_start":           args.stable_id_start,
    }
    prefect_path = f"prefect_params_{args.assembly_name}.json"
    Path(prefect_path).write_text(json.dumps(prefect_params, indent=2))

    # ── Summary ────────────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print(f"  Experiment: {args.accession} ({args.assembly_name})")
    print(f"{'='*65}")
    print(f"  Run script:     {script_path}")
    print(f"  Prefect params: {prefect_path}")
    print(f"  Outdir:         {outdir}")
    print(f"\n  Active stages:")
    print(f"    Ab initio     : Augustus '{augustus_species}' (always)")
    if has_rnaseq:
        print(f"    RNA-seq       : {rnaseq_nf_arg.strip()}")
    if has_genblast:
        print(f"    Protein homol : {uniprot_nf_arg.strip()}")
    if args.refseq_accession:
        print(f"    RefSeq import : {args.refseq_accession}")
    print(f"\n  Layer priorities: {layer_prios_json}")
    print(f"\n  To run:")
    print(f"    tmux new -s {args.assembly_name.lower()}")
    print(f"    salloc --ntasks=1 --mem=4G --time=7-00:00:00 --partition=standard")
    print(f"    bash {script_path}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
