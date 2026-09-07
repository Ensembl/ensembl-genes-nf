# Ensembl stable-ID mapping pipeline

This Nextflow DSL2 pipeline prepares stable-ID updates for an Ensembl core database. It can:

- map stable IDs from a previous live assembly to a newer assembly;
- generate a complete stable-ID reassignment when the current database IDs do not agree with the registry allocation; or
- report that no action is required.

The pipeline generates SQL for review. It does **not** automatically execute the SQL that commits stable-ID changes.

## Requirements

The pipeline is designed for the Ensembl genebuild cluster environment. It expects:

- Nextflow with DSL2 support;
- Python 3 with `pymysql`;
- LiftOn, available as `lifton` by default;
- the `gb1-w` database command for the mapping SQL dry run;
- the `GBS1` and `GBP1` environment variables for the database host and port;
- read access to the target core database and `gb_a_m_test`;
- write access to `gb_a_m_test.annotation_events` when a new mapping session is created;
- access to the Ensembl FTP filesystem trees used by the resolver.

Use the `slurm` profile on the cluster. Input resolution and staging are assigned to the `datamover` queue by this profile.

## Quick start

Run one database and let the pipeline choose the appropriate operation:

```bash
nextflow run main.nf \
    -profile slurm \
    --db_name sus_scrofa_gca000003025v7_core_114_1
```

Resume a previous run after correcting a failure or deleting a published result:

```bash
nextflow run main.nf \
    -profile slurm \
    --db_name sus_scrofa_gca000003025v7_core_114_1 \
    -resume
```

Nextflow reuses a process only when its inputs, script and relevant configuration still match the cached task. Deleting a file under `results/` does not necessarily force its producing process to run again: Nextflow may restore the published file from `work/`. If that process itself must be repeated, use `-resume` after changing its relevant input/script, or remove only that task with Nextflow's cache-management facilities.

## Input modes

Exactly one of `--db_name` or `--samplesheet` is required. They cannot be supplied together.

### Single database

```bash
nextflow run main.nf -profile slurm --db_name DATABASE_NAME
```

Direct file overrides are supported only with `--db_name`:

```bash
nextflow run main.nf \
    -profile slurm \
    --db_name DATABASE_NAME \
    --mode map \
    --ref_fasta /path/reference.fa.gz \
    --ref_gff /path/reference.gff3.gz \
    --target_fasta /path/target.fa.gz \
    --target_gff /path/target.gff3.gz
```

The FASTA and GFF for each assembly are a pair: provide both reference files or neither, and both target files or neither.

### Samplesheet

Example `species.csv`:

```csv
db_name,mode,target_fasta,target_gff,mapping_session_id,ref_fasta,ref_gff
sus_scrofa_gca000003025v7_core_114_1,auto,,,,,
example_species_core,map,/path/target.fa.gz,/path/target.gff3.gz,12345,/path/reference.fa.gz,/path/reference.gff3.gz
another_species_core,reassign,,,,,
```

Run it with:

```bash
nextflow run main.nf -profile slurm --samplesheet species.csv
```

Samplesheet columns:

| Column | Required | Meaning |
|---|---:|---|
| `db_name` | Yes | Target Ensembl core database. Every row must contain it. |
| `mode` | No | `auto`, `map` or `reassign`. A non-empty row value overrides the global `--mode`. |
| `target_fasta` | No | Target assembly FASTA override. Must be accompanied by `target_gff`. |
| `target_gff` | No | Target annotation GFF3 override. Must be accompanied by `target_fasta`. |
| `mapping_session_id` | No | Reuse an existing mapping session instead of creating one. Used only for mapping. |
| `ref_fasta` | No | Reference assembly FASTA override. Must be accompanied by `ref_gff`. |
| `ref_gff` | No | Reference annotation GFF3 override. Must be accompanied by `ref_fasta`. |

Do not combine `--samplesheet` with the command-line file overrides or `--mapping_session_id`. Pipeline-wide parameters such as `--output_dir`, `--rules_config` and `--audit_limit` can still be used.

## Routing modes

### `auto` — default

The resolver uses the database assembly accession and registry information to choose the route:

1. If an earlier version of the same assembly chain is live, run stable-ID mapping from that version.
2. If no live reference exists and current IDs disagree with the registry allocation, generate complete reassignment SQL.
3. If no live reference exists and all IDs agree with the registry allocation, print `NO ACTION REQUIRED` and run neither branch.

### `map`

Force stable-ID mapping. An earlier live assembly version must exist. If it does not, the run fails rather than falling back to reassignment.

If `mapping_session_id` is not supplied, the resolver inserts a new `stable_id_mapping` event in `gb_a_m_test.annotation_events` and uses its ID. A fresh, non-resumed run can therefore create a new mapping session.

### `reassign`

Generate SQL for a complete registry-based reassignment of gene, transcript, translation and exon stable IDs. The reassignment is ordered by each feature table's internal primary key and assigns version `1`.

This route does not require reference or target FASTA/GFF files and does not create a mapping session. If the database already conforms, the generator produces a no-op result rather than changing IDs unnecessarily.

## Input discovery

The resolver reads `assembly.accession` and `species.scientific_name` from the target core database. It obtains stable-ID prefixes and numeric ranges from `gb_a_m_test`.

For a mapping run, omitted files are resolved as follows.

### Reference files

The reference is the newest lower version of the same GCA chain whose genebuild status is `live`. Its files are taken from:

```text
/nfs/ftp/public/ensemblorganisms/<Species_name>/<GCA.version>/genome/unmasked.fa.gz
/nfs/ftp/public/ensemblorganisms/<Species_name>/<GCA.version>/ensembl/geneset/<latest>/genes.gff3.gz
```

### Target files

The resolver first checks the same standard FTP structure. If the target files are unavailable there, it falls back to:

```text
/nfs/ftp/public/databases/ensembl/pre-release/<Species_name>/<GCA.version>/
```

The pre-release directory must contain exactly one `*.dna.softmasked.fa.gz` and exactly one `*.gff3.gz`. Reference files do not use this pre-release fallback.

## Pipeline parameters

### Required choice

| Parameter | Default | Description |
|---|---|---|
| `--db_name` | `null` | Run a single target core database. Mutually exclusive with `--samplesheet`. |
| `--samplesheet` | `null` | CSV containing one or more databases. Mutually exclusive with `--db_name`. |

Exactly one of these parameters must be provided.

### Optional parameters

| Parameter | Default | Applies to | Description |
|---|---|---|---|
| `--mode` | `auto` | Routing | Global mode: `auto`, `map` or `reassign`. A samplesheet row can override it. |
| `--target_fasta` | `null` | Mapping, single database | Target FASTA override. Requires `--target_gff`. |
| `--target_gff` | `null` | Mapping, single database | Target GFF3 override. Requires `--target_fasta`. |
| `--ref_fasta` | `null` | Mapping, single database | Reference FASTA override. Requires `--ref_gff`. |
| `--ref_gff` | `null` | Mapping, single database | Reference GFF3 override. Requires `--ref_fasta`. |
| `--mapping_session_id` | `null` | Mapping, single database | Reuse an existing mapping session. If omitted, a new session is created. |
| `--rules_config` | `assets/stable_id_mapping_rules.json` | Mapping | JSON file containing coordinate and structural matching rules. |
| `--output_dir` | `results` | All routes | Root directory for published pipeline outputs. |
| `--lifton_executable` | `lifton` | Mapping | LiftOn executable or command available to the task. |
| `--lifton_feature_types` | `gene,ncRNA_gene,pseudogene` | Mapping | Comma-separated GFF3 feature types passed to the LiftOn wrapper. |
| `--include_translations` | `true` | Mapping | Include translation stable-ID decisions and SQL. Set to `false` to omit them. |
| `--audit_limit` | `20` | Mapping audit | Maximum number of example rows and some ranked categories printed in the text audit. It does not limit the TSV audit tables. |
| `--replace_events_for_session` | `false` | Mapping SQL | Make the rendered SQL replace existing `stable_id_event` rows for the selected session. It does not replace the `annotation_events` mapping-session row. |
| `--batch_size` | `500` | Mapping and reassignment SQL | Number of rows grouped into generated SQL insert batches. It does not change mapping decisions. |
| `--lifton_threads` | `8` | Currently inactive | Defined in the configuration but not consumed by the active modular workflow. LiftOn currently uses `task.cpus`: 8 in the `slurm` profile and normally 1 without it. |
| `--dry_run_sql` | `false` | Currently inactive | Defined for an older unused module. The active mapping branch always renders executable and dry-run SQL and runs the dry-run SQL. The reassignment branch generates both files but does not execute them. |

Boolean parameters can be changed explicitly, for example:

```bash
nextflow run main.nf \
    -profile slurm \
    --db_name DATABASE_NAME \
    --include_translations false \
    --replace_events_for_session true
```

## Matching rules

The default rules are in `assets/stable_id_mapping_rules.json`. They currently specify:

- minimum coordinate-overlap score: `0.75`;
- structural search window: `100000` bases;
- maximum structural candidates retained: `5`;
- minimum transcript structural score: `0.30`;
- good structural score: `0.45`;
- confident structural score: `0.60`;
- gene transcript-support fraction: `0.60`.

The structural score is weighted mostly toward span containment and query coverage, with smaller contributions from intron similarity, exon agreement, boundary similarity and the LiftOn identity prior. Treat a custom rules file as a substantive mapping-policy change and review it before use.

## How the pipeline works

### Shared resolution and routing

For every database, the pipeline:

1. reads the target assembly and species metadata;
2. obtains the registry stable-ID allocations;
3. finds an earlier live reference assembly when applicable;
4. validates the current database ID populations when reassignment may be needed;
5. resolves files and a mapping session for mapping runs; and
6. routes the database to mapping, reassignment or no action.

Duplicate stable IDs across gene, transcript, translation and exon tables are treated as a fatal error.

### Mapping branch

1. **Stage inputs:** copy or decompress the reference and target FASTA/GFF files into consistently named files.
2. **LiftOn projection:** project the reference annotation onto the target assembly and record reference genes that could not be projected.
3. **Structural matching:** compare projected transcripts with nearby target transcripts. Combine transcript evidence into candidate gene pairs and produce a gene-locus comparison table.
4. **Stable-ID decisions:** attempt structural matching first. If no available structural match is accepted, use coordinate overlap as a fallback. A target feature can be assigned to only one reference feature.
5. **Render SQL:** generate executable SQL and a rollback dry-run version.
6. **Audit:** produce a readable summary plus detailed tables for missing, coordinate-only and new genes.
7. **Database dry run:** execute the rollback SQL through `gb1-w` and save its output.

Mapped genes, transcripts and translations retain their reference stable ID. Unmatched target features receive new IDs from the registry allocation. In the current implementation, target exons receive new assignments rather than participating in retention mapping.

### Reassignment branch

The pipeline checks the complete gene, transcript, translation and exon populations against the registry format and range. When reassignment is required, it generates:

- executable reassignment SQL;
- rollback dry-run SQL; and
- a JSON summary.

The generated reassignment SQL should be reviewed, and the dry-run SQL should be executed and checked before applying the executable SQL.

## Outputs

Outputs are published below `results/<db_name>/` by default.

| Directory | Main contents |
|---|---|
| `inputs/` | Decompressed, consistently named reference and target FASTA/GFF files. |
| `lifton/` | Projected reference GFF3, missing-gene report and LiftOn projection JSON. |
| `matching/` | Transcript pairs, gene pairs, gene-locus comparison and structural-matching JSON. |
| `decisions/` | Stable-ID decisions TSV/JSON and score-evidence TSV. |
| `sql/` | Executable mapping SQL and rollback dry-run SQL. |
| `audit/` | Text audit and detailed missing, coordinate-only and new gene TSV tables. |
| `dry_run_sql/` | Output captured when the mapping dry-run SQL is executed through `gb1-w`. |
| `reassignment/` | Executable reassignment SQL, rollback dry-run SQL and reassignment JSON. Present only for the reassignment route. |

The mapping branch copies staged FASTA and GFF files into `inputs/`; these can be large.

## Reading the mapping audit

The main report is:

```text
results/<db_name>/audit/<db_name>.stable_id_audit.txt
```

### Stable-ID mapping summary

- **Reference:** number of features in the reference annotation.
- **Retained:** reference features whose stable ID was assigned to a target feature.
- **Retained %:** `Retained / Reference`.
- **Missing:** reference features whose stable ID was not retained.
- **Target:** number of features in the target annotation.
- **New:** target features that received newly allocated stable IDs.
- **New %:** `New / Target`.
- **Exons:** currently target-only assignments, so they have no reference-retention percentage.

`Missing` does not necessarily mean that the biological feature disappeared. It means that this run did not accept a one-to-one target assignment for its reference stable ID.

### Version outcomes among retained IDs

This section shows whether the version of each retained stable ID stayed unchanged or was incremented after comparing the relevant old and target feature content. `Unavailable` means that the required comparison could not be performed.

### Mapping evidence among retained IDs

- **Structural:** accepted from projected transcript/exon structure evidence.
- **Coordinate:** structural matching was not used successfully, but the projected and target feature spans passed the coordinate threshold.
- **Other:** a retained mapping whose recorded reason does not belong to the two expected categories.

Coordinate-only mappings deserve more attention because they were accepted using weaker evidence than structural matches.

### Gene annotation classes

- **Reference gene classes:** retention and loss split by GFF feature type and annotation class, such as `protein_coding`, `lncRNA`, `pseudogene` or `IG_V_gene`.
- **Target gene classes:** retained and new target genes split by the same categories.
- **Annotation-class changes:** retained genes whose reference and target classifications differ. These are observations, not automatically errors.

Percentages use the total for the row: reference retention percentages use the reference class count, while target new percentages use the target class count.

### Detailed mapping audit

- **All stable-ID decisions:** counts of `mapped`, `missing` and `new` actions for each feature type.
- **Gene decision audit:** separates structurally mapped genes from coordinate-only mappings, missing reference genes and new target genes.
- **Coordinate-only scores/bands:** distribution of evidence scores for genes rescued by coordinate overlap.
- **Missing gene reasons:** distinguishes genes not projected by LiftOn from genes projected but not matched to an acceptable available target.
- **Missing gene locus status:** compares the best overlapping target gene with the target proposed by structural evidence:
  - `no_locus_candidate`: the projected gene overlaps no target gene on the same sequence and strand;
  - `no_accepted_structure`: an overlapping target exists, but no gene-level structural pair was accepted;
  - `same`: positional and structural evidence point to the same target;
  - `different`: positional and structural evidence point to different targets;
  - `<no locus row>`: normally a gene that LiftOn did not project.
- **Target already claimed:** shows whether a candidate target was assigned to another reference gene. Target assignments are one-to-one.
- **Coordinate-only gene locus status:** explains which structurally unresolved genes were rescued by coordinate overlap.
- **New gene annotation classes:** annotation classes of target genes assigned new stable IDs.
- **New target ID also present in ref GFF:** highlights cases where the target's incoming ID text also occurred in the reference, even though the mapping decision did not retain it automatically.

The counts in these subsections are different classifications of the same decisions and generally must **not** be added together.

### Detailed audit tables

The text report prints only a limited number of examples. The complete rows are in:

- `<db_name>.missing_genes.tsv` — every missing reference gene, its projection/locus evidence and any competing claim;
- `<db_name>.coordinate_mapped_genes.tsv` — genes retained through coordinate fallback and the evidence behind the assignment;
- `<db_name>.new_genes.tsv` — target genes receiving new IDs, their annotation classes and whether they were candidates for old genes.

Use these TSV files for investigation and filtering rather than increasing `audit_limit` excessively.

## Nextflow execution reports

The mapping audit describes the biological/stable-ID results. Nextflow's own reports describe task execution time and resource use.

Use timestamped filenames so repeated runs do not collide:

```bash
run_id=$(date +%Y%m%d_%H%M%S)
mkdir -p results/pipeline_info

nextflow run main.nf \
    -profile slurm \
    --db_name DATABASE_NAME \
    -with-trace results/pipeline_info/trace_${run_id}.txt \
    -with-report results/pipeline_info/report_${run_id}.html \
    -with-timeline results/pipeline_info/timeline_${run_id}.html \
    -with-dag results/pipeline_info/dag_${run_id}.html
```

- **Trace:** one row per task, useful for comparing requested and actual CPU, memory and run time.
- **Report:** interactive summary of execution, resource use and task efficiency.
- **Timeline:** chronological view of when tasks ran and overlapped.
- **DAG:** diagram of the workflow structure.

These are standard Nextflow options, not pipeline parameters; therefore they use a single hyphen.

## Safety and review

Before applying mapping or reassignment SQL:

1. inspect the audit and decision tables;
2. inspect warnings and errors from the dry run;
3. confirm that assignment and matched-row counts agree for every feature type;
4. confirm that the dry-run SQL ends with `ROLLBACK`; and
5. apply executable SQL only after manual review.

The executable SQL files are generated for deliberate manual execution. Do not treat successful pipeline completion alone as approval to modify the database.

