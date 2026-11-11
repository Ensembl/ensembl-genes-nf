#!/usr/bin/env python3
#SBATCH --job-name=dtol-ref-setup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00

"""Prepare DToL long-read references, BAM manifests, and run commands.

The script is deliberately independent of the local workbook and does not
write to a home directory.  It embeds the assembly candidates extracted from
the DToL workbook, checks Ensembl bulk FTP, checks Ensembl pre-release FTP,
falls back to NCBI Datasets, and accepts a candidate only after comparing its
FASTA dictionary with the merged BAM @SQ dictionary.

It can be run directly or submitted with sbatch.  It never submits Nextflow
jobs; it writes an absolute-path run_all_species.sh for an interactive
allocation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


ENSEMBL_BULK = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
ENSEMBL_PRERELEASE = "https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/"
NCBI_DATASETS = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"


# Primary accessions come from DToL_Long_read_Transcriptomics_data_location.
# Older candidates are retained where the earlier BAM-derived map used a
# different accession version; the BAM dictionary decides which is acceptable.
CANDIDATES: dict[str, list[str]] = {
    "Abramis_brama": ["GCA_963993115.1"],
    "Accipiter_gentilis": ["GCA_929443795.2"],
    "Adalia_bipunctata": ["GCA_910592335.1"],
    "Barbatula_barbatula": ["GCA_947034865.1"],
    "Barbus_barbus": ["GCA_936440315.1"],
    "Bombina_variegata": ["GCA_964205305.1"],
    "Calonectris_borealis": ["GCA_964195595.2", "GCA_964195595.1"],
    "Carassius_carassius": ["GCA_963082965.1"],
    "Cetorhinus_maximus": ["GCA_964194155.1"],
    "Coregonus_lavaretus": ["GCA_964263955.1"],
    "Erinaceus_europaeus": ["GCA_950295315.1"],
    "Gasterosteus_aculeatus": ["GCA_964276395.1"],
    "Gobio_gobio": ["GCA_949357685.1"],
    "Gulosus_aristotelis": ["GCA_949628215.1"],
    "Hyperoplus_immaculatus": ["GCA_949357725.1"],
    "Lampetra_fluviatilis": ["GCA_964198595.1"],
    "Larus_michahellis": ["GCA_949357725.1"],
    "Lipophrys_pholis": ["GCA_963383615.1", "GCA_963383615.2"],
    "Lissotriton_vulgaris": ["GCA_964263255.1"],
    "Martes_martes": ["GCA_963455335.1"],
    "Micromesistius_poutassou": ["GCA_964030765.1"],
    "Miniopterus_schreibersii": ["GCA_964146895.1"],
    "Mus_musculus": ["GCA_949316315.1"],
    "Mustelus_asterias": ["GCA_964213995.1"],
    "Myotis_daubentonii": ["GCA_963259705.1", "GCA_963259705.2"],
    "Myotis_mystacinus": ["GCA_964094495.2"],
    "Myotis_nattereri": ["GCA_964212035.1", "GCA_964212035.2"],
    "Myxine_glutinosa": ["GCA_964187855.1"],
    "Numenius_arquata": ["GCA_964106895.1"],
    "Nyctalus_leisleri": ["GCA_964264875.1", "GCA_964264875.2"],
    "Phoxinus_phoxinus": ["GCA_949152265.1"],
    "Pipistrellus_pygmaeus": ["GCA_949987585.1"],
    "Plecotus_auritus": ["GCA_963455305.1"],
    "Pungitius_pungitius": ["GCA_949316345.1"],
    "Rhinolophus_hipposideros": ["GCA_964194185.1"],
    "Rutilus_rutilus": ["GCA_951802725.1"],
    "Scardinius_erythrophthalmus": ["GCA_964197995.1", "GCA_964197995.2"],
    "Scomber_scombrus": ["GCA_963691925.1"],
    "Squalius_cephalus": ["GCA_949319135.1"],
    "Tinca_tinca": [],
    "Triturus_cristatus": ["GCA_964204655.1"],
    "Vespertilio_murinus": ["GCA_963924515.1"],
    "Zootoca_vivipara": ["GCA_963506605.1"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        type=Path,
        default=Path("/hps/nobackup/flicek/ensembl/genebuild/jackt/long-read/dtol_reconstruction"),
        help="Shared output root; must be an absolute path.",
    )
    parser.add_argument(
        "--bam-root",
        type=Path,
        default=Path("/nfs/production/flicek/ensembl/havana/rmerritt/gdm_dtol"),
    )
    parser.add_argument(
        "--pipeline",
        type=Path,
        default=Path("/hps/nobackup/flicek/ensembl/genebuild/jackt/long-read/ensembl-genes-nf"),
    )
    parser.add_argument("--technology", choices=["pacbio", "ont"], default="pacbio")
    parser.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    parser.add_argument("--datasets", type=Path, default=None, help="Existing NCBI datasets executable.")
    parser.add_argument(
        "--no-auto-install-datasets",
        action="store_true",
        help="Fail instead of downloading the NCBI datasets CLI into BASE/bin.",
    )
    parser.add_argument(
        "--allow-unresolved",
        action="store_true",
        help="Return success after preparing the resolved subset; default is strict.",
    )
    return parser.parse_args()


def require_absolute(path: Path, label: str) -> Path:
    path = path.expanduser()
    if not path.is_absolute():
        raise SystemExit(f"{label} must be absolute: {path}")
    return path


def http_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "DToL-long-read-reference-setup/1"})
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read().decode("utf-8", "replace")


def http_exists(url: str) -> bool:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "DToL-long-read-reference-setup/1", "Range": "bytes=0-0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            return response.status in (200, 206)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return False


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "DToL-long-read-reference-setup/1"})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(destination)


def accession_root(accession: str) -> str:
    match = re.fullmatch(r"(GCA|GCF)_(\d{9})\.(\d+)", accession)
    if not match:
        raise ValueError(f"Unsupported assembly accession: {accession}")
    prefix, digits, version = match.groups()
    return f"{ENSEMBL_BULK}{prefix}/{digits[:3]}/{digits[3:6]}/{digits[6:9]}/{version}/ensembl/"


def ensembl_bulk_candidate(species: str, accession: str) -> dict | None:
    root = accession_root(accession)
    try:
        listing = http_text(root)
    except Exception:
        return None

    releases = sorted(set(re.findall(r'href="([0-9]{4}_[0-9]{2})/"', listing)), reverse=True)
    for release in releases:
        release_root = f"{root}{release}/"
        genome_url = ""
        for filename in ("unmasked.fa.bgz", "softmasked.fa.bgz"):
            candidate = f"{release_root}genome/{filename}"
            if http_exists(candidate):
                genome_url = candidate
                break
        if not genome_url:
            continue
        protein_url = f"{release_root}geneset/pep.fa.bgz"
        return {
            "source": "ensembl_bulk",
            "accession": accession,
            "release": release,
            "genome_url": genome_url,
            "protein_url": protein_url if http_exists(protein_url) else "",
        }
    return None


def ensembl_prerelease_candidate(species: str, accession: str) -> dict | None:
    root = f"{ENSEMBL_PRERELEASE}{species}/{accession}/"
    try:
        listing = http_text(root)
    except Exception:
        return None

    files = re.findall(r'href="([^"]+)"', listing)
    genome_names = [
        name for name in files
        if re.search(r"\.dna\.(unmasked|softmasked)\.fa(?:\.gz)?$", name)
    ]
    genome_names.sort(key=lambda name: ("unmasked" not in name, name))
    if not genome_names:
        return None

    protein_names = [name for name in files if re.search(r"\.pep(?:\.all)?\.fa(?:\.gz)?$", name)]
    protein_names.sort()
    genome_url = urllib.parse.urljoin(root, genome_names[0])
    protein_url = urllib.parse.urljoin(root, protein_names[0]) if protein_names else ""
    return {
        "source": "ensembl_prerelease",
        "accession": accession,
        "release": "pre-release",
        "genome_url": genome_url,
        "protein_url": protein_url,
    }


def ensure_datasets(args: argparse.Namespace, base: Path) -> Path:
    if args.datasets:
        executable = require_absolute(args.datasets, "--datasets")
        if not executable.is_file():
            raise SystemExit(f"NCBI datasets executable not found: {executable}")
        return executable

    found = shutil.which("datasets")
    if found:
        return Path(found).resolve()

    if args.no_auto_install_datasets:
        raise SystemExit("NCBI datasets is required for the fallback but was not found in PATH")

    executable = base / "bin" / "datasets"
    if not executable.exists():
        print(f"Installing NCBI datasets under {executable}")
        download(NCBI_DATASETS, executable)
        executable.chmod(executable.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return executable


def ncbi_candidate(species: str, accession: str, datasets: Path, base: Path) -> dict | None:
    root = base / "ncbi" / species / accession
    package_zip = root / "ncbi_dataset.zip"
    package_dir = root / "package"
    root.mkdir(parents=True, exist_ok=True)

    if not package_dir.exists():
        if not package_zip.exists():
            command = [
                str(datasets),
                "download",
                "genome",
                "accession",
                accession,
                "--include",
                "genome,protein,gff3,gtf,seq-report",
                "--filename",
                str(package_zip),
                "--no-progressbar",
            ]
            print("NCBI:", " ".join(command))
            try:
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError:
                return None
        try:
            with zipfile.ZipFile(package_zip) as archive:
                archive.extractall(package_dir)
        except (OSError, zipfile.BadZipFile):
            return None

    fasta_candidates = [
        path for path in package_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".fna", ".fa", ".fasta"}
    ]
    if not fasta_candidates:
        return None
    fasta_candidates.sort(key=lambda path: ("genomic" not in path.name, -path.stat().st_size))

    protein_candidates = [
        path for path in package_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".faa", ".fa", ".fasta"} and "protein" in path.name.lower()
    ]
    protein_candidates.sort(key=lambda path: -path.stat().st_size)

    return {
        "source": "ncbi_datasets",
        "accession": accession,
        "release": "NCBI",
        "genome_path": str(fasta_candidates[0]),
        "protein_path": str(protein_candidates[0]) if protein_candidates else "",
    }


def materialize(source: str, destination: Path, cache: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        return

    source_path: Path | None = None
    if source.startswith("http://") or source.startswith("https://"):
        parsed = urllib.parse.urlparse(source)
        source_path = cache / Path(parsed.path).name
        if not source_path.exists() or source_path.stat().st_size == 0:
            download(source, source_path)
    else:
        source_path = Path(source)

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    if source_path.name.endswith((".gz", ".bgz")):
        with gzip.open(source_path, "rb") as source_handle, partial.open("wb") as destination_handle:
            shutil.copyfileobj(source_handle, destination_handle, length=1024 * 1024)
    else:
        shutil.copyfile(source_path, partial)
    partial.replace(destination)


def bam_dictionary(bam: Path) -> dict[str, int]:
    header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
    result: dict[str, int] = {}
    for line in header.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = {}
        for field in line.rstrip().split("\t")[1:]:
            if ":" in field:
                key, value = field.split(":", 1)
                fields[key] = value
        if "SN" in fields and "LN" in fields:
            result[fields["SN"]] = int(fields["LN"])
    return result


def fasta_dictionary(fasta: Path) -> dict[str, int]:
    fai = Path(str(fasta) + ".fai")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    result: dict[str, int] = {}
    with fai.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            result[fields[0]] = int(fields[1])
    return result


def compare_dictionaries(bam_sq: dict[str, int], ref_sq: dict[str, int]) -> tuple[str, str]:
    missing = sorted(set(bam_sq) - set(ref_sq))
    mismatched = sorted(name for name in bam_sq if name in ref_sq and bam_sq[name] != ref_sq[name])
    if missing or mismatched:
        return "mismatch", f"missing={','.join(missing)} length_mismatch={','.join(mismatched)}"
    if list(bam_sq.items()) == list(ref_sq.items()):
        return "exact", "all BAM contigs and lengths match"
    if set(bam_sq) == set(ref_sq):
        return "mismatch", "BAM and reference contigs have the same names and lengths but a different order"
    return "compatible", f"BAM contigs match; reference has {len(set(ref_sq) - set(bam_sq))} extra contigs"


def ordered_reference_mapping(
    bam_sq: dict[str, int], ref_sq: dict[str, int]
) -> tuple[dict[str, str], str] | None:
    """Return a safe reference-name -> BAM-name mapping, if one exists.

    DToL BAMs can retain the accession names from the assembly used for
    alignment, while Ensembl FASTAs commonly use chromosome names (for
    example ``1`` instead of ``OZ023076.1``).  We only infer a rename when
    the BAM dictionary can be recovered as an ordered subsequence of the
    reference dictionary with identical lengths.  This preserves sequence
    order and prevents an arbitrary length-only match from changing the
    assembly coordinates.
    """

    bam_names = list(bam_sq)
    ref_names = list(ref_sq)

    if list(bam_sq.items()) == list(ref_sq.items()):
        return None

    # A reference may contain unplaced/alternate contigs in addition to the
    # primary contigs retained by the BAM.  If the BAM names are already
    # present, dropping those extras is safe only when the retained order is
    # identical to the BAM order and all lengths agree.
    if all(name in ref_sq and ref_sq[name] == bam_sq[name] for name in bam_names):
        retained = [name for name in ref_names if name in bam_sq]
        if retained == bam_names:
            return {name: name for name in bam_names}, "dropped reference-only contigs"
        if set(retained) == set(bam_names):
            return None

    # Find a unique ordered subsequence of reference contigs whose lengths
    # equal the BAM lengths.  Counting is capped at two: we only need to know
    # whether the solution is unique or ambiguous.
    n_ref = len(ref_names)
    n_bam = len(bam_names)
    counts = [[0] * (n_bam + 1) for _ in range(n_ref + 1)]
    for i in range(n_ref + 1):
        counts[i][n_bam] = 1
    for i in range(n_ref - 1, -1, -1):
        for j in range(n_bam - 1, -1, -1):
            count = counts[i + 1][j]
            if ref_sq[ref_names[i]] == bam_sq[bam_names[j]]:
                count += counts[i + 1][j + 1]
            counts[i][j] = min(count, 2)

    if counts[0][0] != 1:
        return None

    mapping: dict[str, str] = {}
    i = 0
    for bam_name in bam_names:
        while i < n_ref:
            ref_name = ref_names[i]
            i += 1
            if ref_sq[ref_name] != bam_sq[bam_name]:
                continue
            if counts[i][len(mapping) + 1] == 0:
                continue
            mapping[ref_name] = bam_name
            break
        else:
            return None
    return mapping, "renamed headers using unique ordered contig lengths"


def normalize_reference_to_bam(
    reference: Path,
    bam_sq: dict[str, int],
    output: Path,
    mapping_path: Path,
) -> tuple[bool, str]:
    """Rewrite FASTA headers to the BAM dictionary when the mapping is safe."""

    ref_sq = fasta_dictionary(reference)
    mapping_result = ordered_reference_mapping(bam_sq, ref_sq)
    if mapping_result is None:
        return False, "no unique ordered contig-length mapping"
    mapping, method = mapping_result

    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    seen: set[str] = set()
    current_keep = False
    with reference.open() as source, partial.open("w") as destination:
        for line in source:
            if line.startswith(">"):
                source_name = line[1:].split(None, 1)[0]
                target_name = mapping.get(source_name)
                current_keep = target_name is not None
                if current_keep:
                    if target_name in seen:
                        partial.unlink(missing_ok=True)
                        return False, f"duplicate mapped BAM contig {target_name}"
                    seen.add(target_name)
                    destination.write(f">{target_name}\n")
                continue
            if current_keep:
                destination.write(line)
    partial.replace(output)
    Path(str(output) + ".fai").unlink(missing_ok=True)

    normalized_sq = fasta_dictionary(output)
    match_kind, detail = compare_dictionaries(bam_sq, normalized_sq)
    if match_kind != "exact":
        return False, f"normalized FASTA failed exact validation: {detail}"

    with mapping_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_name", "bam_name", "length"])
        for source_name, bam_name in mapping.items():
            writer.writerow([source_name, bam_name, bam_sq[bam_name]])
    return True, method


def discover_bam(bam_root: Path, species: str) -> Path | None:
    directory = bam_root / species / "reads" / "all_mapped_reads" / "merged"
    bams = sorted(directory.glob("*.bam"))
    if len(bams) == 1:
        return bams[0]
    if not bams:
        return None
    raise RuntimeError(f"{species}: expected one merged BAM, found {len(bams)} in {directory}")


def stage_bam(source: Path, destination: Path, threads: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or destination.exists():
        destination.unlink()
    destination.symlink_to(source)

    source_indexes = [
        Path(str(source) + ".csi"),
        source.with_suffix(".csi"),
        Path(str(source) + ".bai"),
        source.with_suffix(".bai"),
    ]
    source_index = next((path for path in source_indexes if path.exists()), None)
    destination_suffix = ".csi" if source_index is None else source_index.suffix
    destination_index = Path(str(destination) + destination_suffix)
    if destination_index.is_symlink() or destination_index.exists():
        destination_index.unlink()
    if source_index:
        destination_index.symlink_to(source_index)
    else:
        subprocess.run(
            ["samtools", "index", "-c", "-@", str(threads), "-o", str(destination_index), str(source)],
            check=True,
        )
    subprocess.run(["samtools", "quickcheck", "-v", str(destination)], check=True)


def write_manifest(
    path: Path,
    species: str,
    accession: str,
    source_bam: Path,
    staged_bam: Path,
    technology: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "id": species,
        "cohort_id": f"dtol_{species}",
        "run_id": species,
        "technology": technology,
        "minimap2_preset": "splice:hq" if technology == "pacbio" else "splice",
        "isoquant_data_type": "pacbio_ccs" if technology == "pacbio" else "nanopore",
        "transgram_data_type": "auto",
        "assembly_id": accession,
        "read_set_id": "all_mapped_reads",
        "sample_id": species,
        "library_id": "dtol",
        "read_representation": "bam",
        "source_uri": str(source_bam),
        "bam": str(staged_bam),
    }
    fields = list(row)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\\''") + "'"


def write_run_script(base: Path, pipeline: Path, resolutions: list[dict]) -> Path:
    path = base / "run_all_species.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"export BASE={shell_quote(str(base))}",
        f"export PIPELINE={shell_quote(str(pipeline))}",
        "export NXF_HOME=\"$BASE/nxf-home\"",
        "export NXF_WORK=\"$BASE/work\"",
        "export NXF_SINGULARITY_CACHEDIR=\"$BASE/singularity-cache\"",
        "export SINGULARITY_CACHEDIR=\"$BASE/singularity-cache\"",
        "export TMPDIR=\"$BASE/tmp\"",
        "",
        "mkdir -p \"$NXF_HOME\" \"$NXF_WORK\" \"$NXF_SINGULARITY_CACHEDIR\" \"$TMPDIR\"",
        "",
    ]

    for item in resolutions:
        species = item["species"]
        outdir = base / "results" / species
        workdir = base / "work" / species
        manifest = base / "manifests" / f"{species}.tsv"
        reference = base / "references" / species / "reference.fa"
        proteins = base / "references" / species / "proteins.faa"
        lines.extend(
            [
                f"mkdir -p {shell_quote(str(outdir))}",
                "",
                "nextflow run \"$PIPELINE/pipelines/long-read/main.nf\" \\",
                "    -c \"$PIPELINE/config/singularity.config\" \\",
                "    -process.executor local \\",
                f"    -work-dir {shell_quote(str(workdir))} \\",
                "    --mode bam \\",
                f"    --manifest {shell_quote(str(manifest))} \\",
                f"    --reference_fasta {shell_quote(str(reference))} \\",
                f"    --cohort_id {shell_quote('dtol_' + species)} \\",
                "    --model_backends all \\",
                "    --backend_scopes run \\",
                "    --isoquant_scope run \\",
                "    --bambu_scope run \\",
                "    --shard_max_concurrent 1 \\",
                "    --evidence_mode enabled \\",
                "    --orf_mode enabled \\",
            ]
        )
        if proteins.exists() and proteins.stat().st_size > 0:
            lines.extend([
                "    --diamond_mode enabled \\",
                f"    --diamond_reference_proteins {shell_quote(str(proteins))} \\",
            ])
        else:
            lines.append("    --diamond_mode disabled \\")
        lines.extend([
            f"    --outdir {shell_quote(str(outdir))} \\",
            f"    --tracedir {shell_quote(str(outdir / 'pipeline_info'))} \\",
            f"    -with-report {shell_quote(str(outdir / 'pipeline_report.html'))} \\",
            f"    -with-timeline {shell_quote(str(outdir / 'pipeline_timeline.html'))} \\",
            "    -resume",
            "",
        ])

    path.write_text("\n".join(lines))
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return path


def main() -> int:
    args = parse_args()
    base = require_absolute(args.base, "--base")
    bam_root = require_absolute(args.bam_root, "--bam-root")
    pipeline = require_absolute(args.pipeline, "--pipeline")

    for executable in ("samtools", "bgzip"):
        if not shutil.which(executable):
            raise SystemExit(f"Required executable not found in PATH: {executable}")
    if not (pipeline / "pipelines/long-read/main.nf").is_file():
        raise SystemExit(f"Long-read pipeline not found under {pipeline}")

    for directory in ("metadata", "references", "bams", "manifests", "results", "work", "tmp", "download_cache"):
        (base / directory).mkdir(parents=True, exist_ok=True)

    species_with_bams = []
    for child in sorted(bam_root.iterdir()):
        if child.is_dir() and (child / "reads" / "all_mapped_reads" / "merged").is_dir():
            species_with_bams.append(child.name)

    if not species_with_bams:
        raise SystemExit(f"No DToL merged BAM directories found below {bam_root}")

    datasets = None

    resolutions: list[dict] = []
    failures: list[dict] = []

    for species in species_with_bams:
        print(f"\n=== {species} ===")
        source_bam = discover_bam(bam_root, species)
        assert source_bam is not None
        bam_sq = bam_dictionary(source_bam)
        candidates = CANDIDATES.get(species, [])
        if not candidates:
            failures.append({"species": species, "reason": "no assembly accession encoded"})
            print("FAIL: no assembly accession encoded")
            continue

        matches: list[dict] = []
        for accession in candidates:
            sources = []
            bulk = ensembl_bulk_candidate(species, accession)
            if bulk:
                sources.append(bulk)
            prerelease = ensembl_prerelease_candidate(species, accession)
            if prerelease:
                sources.append(prerelease)

            # Only install/find NCBI Datasets when neither Ensembl location
            # contains this exact accession.  This keeps ordinary Ensembl-only
            # preparation lightweight while making the fallback self-contained.
            if not sources:
                if datasets is None:
                    datasets = ensure_datasets(args, base)
                ncbi = ncbi_candidate(species, accession, datasets, base)
                if ncbi:
                    sources.append(ncbi)

            for candidate in sources:
                candidate_cache = base / "candidate_references" / species / accession / candidate["source"]
                candidate_cache.mkdir(parents=True, exist_ok=True)
                candidate_fasta = candidate_cache / "reference.fa"

                genome_source = candidate.get("genome_path") or candidate.get("genome_url")
                try:
                    materialize(genome_source, candidate_fasta, base / "download_cache")
                    ref_sq = fasta_dictionary(candidate_fasta)
                    match_kind, detail = compare_dictionaries(bam_sq, ref_sq)
                    header_map = None
                    if match_kind != "exact":
                        normalized_fasta = candidate_cache / "reference.bam-named.fa"
                        header_map = candidate_cache / "reference_header_map.tsv"
                        normalized, normalization_detail = normalize_reference_to_bam(
                            candidate_fasta,
                            bam_sq,
                            normalized_fasta,
                            header_map,
                        )
                        if normalized:
                            candidate_fasta = normalized_fasta
                            match_kind = "normalized_exact"
                            detail = normalization_detail
                        else:
                            detail = f"{detail}; {normalization_detail}"
                except (OSError, subprocess.CalledProcessError, urllib.error.URLError) as error:
                    print(f"{accession} {candidate['source']}: validation error: {error}")
                    continue

                print(f"{accession} {candidate['source']}: {match_kind} ({detail})")
                if match_kind in {"exact", "normalized_exact"}:
                    candidate["match_kind"] = match_kind
                    candidate["match_detail"] = detail
                    candidate["candidate_fasta"] = candidate_fasta
                    candidate["header_map"] = header_map
                    matches.append(candidate)

        if not matches:
            failures.append({"species": species, "reason": "no candidate reference matched BAM @SQ dictionary"})
            print("FAIL: no candidate reference matched BAM @SQ dictionary")
            continue

        score = {"exact": 2, "normalized_exact": 1}
        best_score = max(score[item["match_kind"]] for item in matches)
        best = [item for item in matches if score[item["match_kind"]] == best_score]
        accessions = sorted({item["accession"] for item in best})
        if len(accessions) > 1:
            failures.append({"species": species, "reason": f"ambiguous matching accessions: {','.join(accessions)}"})
            print(f"FAIL: ambiguous matching accessions: {', '.join(accessions)}")
            continue
        selected = best[0]

        ref_dir = base / "references" / species
        ref_dir.mkdir(parents=True, exist_ok=True)
        reference = ref_dir / "reference.fa"
        shutil.copyfile(selected["candidate_fasta"], reference)
        ref_sq = fasta_dictionary(reference)
        final_match_kind, final_detail = compare_dictionaries(bam_sq, ref_sq)
        if final_match_kind != "exact":
            failures.append({"species": species, "reason": f"selected reference failed final BAM dictionary validation: {final_detail}"})
            print(f"FAIL: selected reference failed final BAM dictionary validation: {final_detail}")
            continue

        header_map = selected.get("header_map")
        final_header_map = ref_dir / "reference_header_map.tsv"
        if header_map and Path(header_map).is_file():
            shutil.copyfile(header_map, final_header_map)

        protein_source = selected.get("protein_path") or selected.get("protein_url")
        protein = ref_dir / "proteins.faa"
        if protein_source:
            try:
                materialize(protein_source, protein, base / "download_cache")
            except (OSError, urllib.error.URLError) as error:
                print(f"Protein download failed for {species}: {error}")

        staged_bam = base / "bams" / f"{species}.bam"
        stage_bam(source_bam, staged_bam, args.threads)
        manifest = base / "manifests" / f"{species}.tsv"
        write_manifest(manifest, species, selected["accession"], source_bam, staged_bam, args.technology)

        metadata = {
            "species": species,
            "assembly_accession": selected["accession"],
            "source": selected["source"],
            "release": selected.get("release", ""),
            "genome_url": selected.get("genome_url", ""),
            "protein_url": selected.get("protein_url", ""),
            "bam": str(source_bam),
            "match_kind": selected["match_kind"],
            "match_detail": selected["match_detail"],
            "header_map": str(final_header_map) if final_header_map.exists() else "",
        }
        (ref_dir / "reference_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        resolutions.append({**metadata, "reference": str(reference), "protein": str(protein) if protein.exists() else ""})
        print(f"SELECTED {selected['accession']} via {selected['source']} ({selected['match_kind']})")

    resolution_path = base / "metadata" / "reference_resolution.tsv"
    with resolution_path.open("w", newline="") as handle:
        fields = ["species", "assembly_accession", "source", "release", "match_kind", "match_detail", "header_map", "reference", "protein"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in resolutions:
            writer.writerow({field: item.get(field, "") for field in fields})

    failure_path = base / "metadata" / "reference_failures.tsv"
    with failure_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["species", "reason"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(failures)

    run_script = write_run_script(base, pipeline, resolutions)
    print(f"\nResolved species: {len(resolutions)}")
    print(f"Unresolved species: {len(failures)}")
    print(f"Resolution report: {resolution_path}")
    print(f"Failure report: {failure_path}")
    print(f"Nextflow run script: {run_script}")

    if failures and not args.allow_unresolved:
        print("Strict mode: unresolved species remain; no full-species run is authorised.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
