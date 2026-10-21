# Track-file preparation

This DSL2 pipeline treats the supplied TSV as the authoritative entity
inventory. It does not scan production directories or infer run/sample/merged
relationships. The per-row `track_types` value is used unless the pipeline
parameter supplies an explicit override; `all` means coverage, gene-model,
and splice-junction tracks. Junction derivation from BAM is opt-in with
`--derive_junctions true`.

Run-level manifests can be discovered with the standalone generator. `main`
and `anno` default to all three tracks; `gdm` defaults to coverage and
splice-junction because GDM input is handled at run level and tissue GTFs are
out of scope.

```bash
python3 pipelines/track_file_preparation/bin/generate_manifest.py \
  --mode main \
  --base-dir /path/to/GCA_052040795.2 \
  --output-manifest /path/to/manifests/main.tsv

python3 pipelines/track_file_preparation/bin/generate_manifest.py \
  --mode gdm \
  --base-dir /lts/production/fergal/genebuild/gdm_data \
  --gdm-duplicate-policy newest \
  --output-manifest /path/to/manifests/gdm_runs.tsv
```

`gdm` searches for accession-named STAR BAMs below `alignment/` and matches
their adjacent SJ files. Its base can be a single species/taxon directory or
the whole `gdm_data` tree; GCA accessions are read from canonical GCA path
components or species directory names such as `sus_scrofa_gca000003025v6`. For
a multi-species root, omit `--assembly-release` so each row defaults to its own
GCA release.
GDM manifests contain run entities only. Tissue GTFs are intentionally not
processed by this workflow.

If a run accession appears more than once within the same assembly, the default
`--gdm-duplicate-policy newest` keeps the uniquely newest BAM; the selected row
records that choice in the audit warning. The same run accession in different
assemblies produces separate rows. Tied mtimes remain unresolved. When both
`.bai` and `.csi` exist, the manifest uses `.bai`.

When any manifest path is under `/lts`, supply a persistent HPS/NFS destination
with `--lts_stage_dir`. The pipeline copies each entity's LTS inputs there at
the start, on the Slurm `datamover` partition, preserving source timestamps.
For example:

```bash
--lts_stage_dir /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/staged_gdm
```

The stage directory keeps copies named by assembly, entity, and input type.
Non-LTS paths continue to be used in place. An index older than its BAM is
discarded and regenerated from the staged BAM before BAM-dependent tracks run.

Before BAM-based tracks run, the pipeline compares each BAM sequence dictionary
with the supplied `chrom.sizes`. If names and lengths match, it uses a sorted
copy of the supplied file. If names differ but the complete set of sequence
lengths matches, it writes a BAM-named `chrom.sizes` for that entity; this
supports alignment aliases such as `1` where the assembly file uses accession
IDs. A length mismatch stops the run before track generation. The resolution
and chromosome count are recorded under `${outdir}/reference/`.

The generator writes an audit, exceptions table, and JSON summary beside the
manifest. Strict mode is the default: any unresolved candidate prevents the
manifest from being published. `--dry-run` writes diagnostics without the
manifest.

Run the wiring test with:

```bash
nextflow run pipelines/track_file_preparation/main.nf -stub-run -profile test \
  --input_manifest pipelines/track_file_preparation/test/entity_manifest.tsv \
  --chrom_sizes pipelines/track_file_preparation/test/test.chrom.sizes \
  --assembly_release test-release --outdir /tmp/track-file-preparation-stub
```

The manifest must contain the required columns from
`TRACK_FILE_PREPARATION_SPEC.md`. Relative paths are resolved against the
manifest directory or `--manifest_base_dir`. Output records retain the
original entity ID while filenames use a collision-checked sanitized ID.

For BAM-dependent tracks, the pipeline uses a manifest-supplied BAM index or
an adjacent `.bai`/`.csi` only when its modification time is at least as new
as the BAM. If no current index exists, it creates a `.bai` with `samtools index`
before running coverage or BAM-derived junction extraction. Index
creation is a declared workflow stage, not a prerequisite that the manifest
producer must perform.
