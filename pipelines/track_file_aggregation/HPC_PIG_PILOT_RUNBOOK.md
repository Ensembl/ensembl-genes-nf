# HPC pig pilot aggregation runbook

This records the direct run-to-tissue aggregation commands for the 200-run
`GCA_000003025.6` pilot. It assumes run-level preparation has completed and
that `run_level/track_manifest.tsv` contains complete coverage and
splice-junction records. Gene-model BigBeds are not aggregated.

## Create the explicit aggregation manifest

Run on `hl-codon-dm-01`:

```bash
python3 - <<'PY'
import csv
import json
from pathlib import Path

metadata_path = Path(
    "/hps/nobackup/flicek/ensembl/genebuild/jackt/"
    "expression_tracks/pig-pilot/"
    "sus_scrofa_gca000003025v6_representative_200.metadata.tsv"
)

output_path = Path(
    "/hps/nobackup/flicek/ensembl/genebuild/jackt/"
    "expression_tracks/pig-pilot/run_to_tissue_aggregation.tsv"
)

priority = [
    "brain_cns", "muscle", "liver", "kidney", "heart", "lung_airway",
    "blood_immune", "gut_intestine", "reproductive", "embryo_fetus",
    "skin_adipose", "cell_line_culture",
]

def choose_tissue(raw):
    tags = set(raw.split(","))
    for tissue in priority:
        if tissue in tags:
            return tissue
    return "uncharacterized"

with metadata_path.open() as source, output_path.open("w", newline="") as destination:
    reader = csv.DictReader(source, delimiter="\t")
    fields = [
        "aggregate_id", "aggregate_type", "source_entity_id",
        "source_entity_type", "gca_accession", "aggregate_label",
        "aggregation_method", "metadata_json",
    ]
    writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for row in reader:
        tissue = choose_tissue(row["tissue_tags"])
        writer.writerow({
            "aggregate_id": f"tissue_{tissue}",
            "aggregate_type": "tissue",
            "source_entity_id": row["run_accession"],
            "source_entity_type": "run",
            "gca_accession": "GCA_000003025.6",
            "aggregate_label": tissue,
            "aggregation_method": "raw",
            "metadata_json": json.dumps({
                "tissue": tissue,
                "metadata_source": "ENA sample title/description heuristic",
            }, sort_keys=True),
        })
PY
```

## Run tissue aggregation

```bash
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/ensembl-genes-nf

nextflow run \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/ensembl-genes-nf/pipelines/track_file_aggregation/main.nf \
  -profile slurm \
  -resume \
  -name pig_pilot_tissue_aggregates_200 \
  --track_manifest \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/run_level/track_manifest.tsv \
  --aggregation_manifest \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/run_to_tissue_aggregation.tsv \
  --manifest_base_dir \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/run_level \
  --chrom_sizes \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/reference/Sscrofa11.1/GCA_000003025.6.ensembl.chrom.sizes \
  --track_types all \
  --outdir \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/tissue_aggregates \
  -work-dir \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/work/tissue_aggregates \
  -with-report \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/tissue_aggregates/nextflow_report.html \
  -with-trace \
  /hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/tissue_aggregates/tissue_aggregates_trace.txt
```
