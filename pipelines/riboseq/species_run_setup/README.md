# RDP production run setup

Use the local RiboSeq Data Portal release SQLite database as the run inventory.
The exporter opens it read-only and records the SHA-256 of the release, the
reference manifest, the exporter and any optional review or selection input.
It publishes a complete bundle to a new directory after the inputs and quotas
validate.

## Export a release

From the `ensembl-genes-nf` checkout with Python 3.9 or later:

```bash
RDP=/Users/jackt/projects/RiboSeqOrg-DataPortal/build/metadata_rebuild/release_2026-11/output
BUNDLE="$PWD/results/rdp-production/release_2026-11-v4"
python3 pipelines/riboseq/bin/prepare_rdp_processing.py \
  --db "$RDP/db.sqlite3" \
  --review-queue "$RDP/review_queue.csv" \
  --pilot-plan pipelines/riboseq/species_run_setup/rdp_pilot_plan.json \
  --out-dir "$BUNDLE"
```

Use a fresh, versioned destination when regenerating an export. The SQLite
release and its optional review queue stay local. Stage the generated bundle
on the execution filesystem when launching on Codon.

| Output | Use |
| --- | --- |
| `run_inventory.tsv` | One row per portal accession with curated library type, source metadata, reference ID, hold reasons and raw-input status. |
| `species_summary.tsv` | Counts by portal scientific name and TaxID, including pending reference mappings. |
| `references.tsv` and `organism_setup/*.json` | Pinned assembly, geneset and current file URLs with setup-only Nextflow parameters. |
| `sample_sheets/*_samplesheet.csv` | Full candidate list per configured species. |
| `controls/*_samplesheet.csv` and `pilot/*_samplesheet.csv` | Named-control and draft pilot sheets. The 12 controls are included in the 225 human, 225 mouse and 50 zebrafish quotas. |
| `pilot/selection.tsv` and `provenance.json` | Named control intent, selection strata, seed, source hashes and release totals. |
| `processing.json` | Per-sample getRPF settings and exact single-end ENA acquisition. |

Pass the `*_samplesheet.csv` files directly to `--sample_sheet`. They follow
the checked-in species sheets with exactly two columns, `Run,study_accession`,
and CRLF line endings. `Run` is an archive run accession and
`study_accession` is its BioProject accession.

The baseline uses the portal's curated `Ribo-Seq` label, `SINGLE` layout and
a BioProject. It holds unresolved review-queue rows, other assays and layouts,
and taxa without a configured mapping. Reference joins require both a TaxID
listed in the manifest and a matching scientific name. Review unmapped taxa
and strain assignments before expanding the reference catalog. Portal
processing status is historical evidence. All raw inputs start as unresolved.

The pilot selection cycles through studies and spreads each round across
inhibitor, nuclease, raw read length and depth bins. Missing measurements have
their own bins. This is a study-balanced draft for evaluating coverage across
protocols. Keep the named controls visible when calculating routine-stratum
success rates.

## Build references and run the controls

The v4 export contains the GRCz13ab reference choice. Stage it on Codon at
`/homes/jackt/rdp-production/inputs/release_2026-11-v4`. Set `OUT` to a
writable execution directory for reference and control outputs.

For mouse, reuse the matching GRCm39 / Ensembl `2026_04` bundle with
`/homes/jackt/rdp-production/reference-configs/release_2026-11-v2/mus_musculus.config`.
Its paths point to the verified, flat October reference directory. Skip the
setup command and pass this file with `-c` for the mouse controls.

Build the pinned GRCh38.p14 / Ensembl `2026_04` human and GRCz13ab
pre-release zebrafish reference with their setup JSONs. The zebrafish files
are the soft-masked FASTA and GTF under the Ensembl pre-release directory.
Use the pinned workflow checkout and a persistent Nextflow work directory.
When launching from the repository root, pass the Ribo-seq pipeline config as
a global `-c` option before `run`; this is what enables the Singularity profile
for the `gffread`, Bowtie, RiboMetric and version-collation containers.

```bash
export NXF_VER=26.04.6
BUNDLE=/homes/jackt/rdp-production/inputs/release_2026-11-v4
SPECIES=danio_rerio
REF="$OUT/references/$SPECIES"
/hps/software/spack/opt/spack/linux-rhel9-zen5/none-none/nextflow-26.04.0-5if3abirupmf5qdo2jpsqte52235gjpy/bin/nextflow \
  -c /hps/software/users/ensembl/genebuild/jackt/ensembl-genes-nf/pipelines/riboseq/nextflow.config \
  run pipelines/riboseq/main.nf -profile slurm \
  -params-file "$BUNDLE/organism_setup/$SPECIES.json" \
  -work-dir "$OUT/work/references/$SPECIES" --outdir "$REF"

/hps/software/spack/opt/spack/linux-rhel9-zen5/none-none/nextflow-26.04.0-5if3abirupmf5qdo2jpsqte52235gjpy/bin/nextflow \
  -c /hps/software/users/ensembl/genebuild/jackt/ensembl-genes-nf/pipelines/riboseq/nextflow.config \
  -c "$REF/pipeline_info/riboseq_params.config" \
  run pipelines/riboseq/main.nf -profile slurm \
  -params-file "$BUNDLE/processing.json" \
  -work-dir "$OUT/work/controls/$SPECIES" \
  --sample_sheet "$BUNDLE/controls/${SPECIES}_samplesheet.csv" \
  --outdir "$OUT/controls/$SPECIES"
```

Capture reference checksums and the resolved pipeline config with each
attempt. Exact acquisition uses the ENA file report for the requested run,
checks that the archive reports one single-end FASTQ, and verifies its MD5
before downstream processing. The successful source record, including a
local SHA-256, is published under `pipeline_info/fastq_sources/`. Missing
ENA FASTQs, ambiguous files and unsupported layouts are acquisition failures
with the accession and reason in the task log. The exact-acquisition mode
currently requires an ENA-hosted FASTQ. Keep accessions that need another
acquisition route in the review ledger.

Look up a control's source without downloading its library:

```bash
python3 pipelines/riboseq/bin/fetch_ena_run.py SRR11798519 --metadata-only
```

Use the control results to close the acquisition, trimming, contaminant,
alignment-count and P-site checks before freezing the pilot. The generated
processing parameters use per-sample mode.

## Reference snapshot

`species_manifest.tsv` preserves the earlier `FINAL_LIST.csv` discovery
snapshot and its pinned reference choices. The checked-in `sample_sheets/`
and `configs/` next to it record that historical setup. The exporter carries
forward the assembly and geneset, translating retired species-based FTP paths
to Ensembl's accession-based tree. It uses the unmasked genome and
`genes.gtf.gz` annotations. [Ensembl FTP migration](https://www.ensembl.info/2026/06/)
