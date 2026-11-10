# Protein database pipeline

This pipeline builds reproducible DIAMOND protein databases from taxon-scoped
UniProtKB or UniProtKB/Swiss-Prot FASTA records. It is deliberately separate
from the long-read pipeline: the long-read workflow consumes the published
database.dmnd or database_proteins.faa path and does not fetch external
reference data during a production run.

UniProtKB and Swiss-Prot are selected through the build manifest:

- uniprot retrieves the broader UniProtKB set, including reviewed Swiss-Prot
  and unreviewed TrEMBL entries.
- swissprot adds reviewed:true, restricting the query to the reviewed
  UniProtKB/Swiss-Prot subset.

The scope column controls the taxonomy query. species uses an exact NCBI
taxonomy ID; clade uses the UniProt lineage query field and includes descendant
taxa. Each build records the exact query, API endpoint, reported UniProt
release, retrieval time, sequence count, FASTA checksum, DIAMOND version, and
database checksums.

## Build manifest

The manifest is tab-separated and has these columns:

| Column | Required | Meaning |
| --- | --- | --- |
| build_id | yes | Unique output-safe identifier; becomes the output directory name |
| source | yes | uniprot or swissprot |
| scope | yes | species or clade |
| taxon_id | yes | NCBI taxonomy ID |
| taxon_name | no | Human-readable label stored in provenance |
| query_extra | no | Optional additional UniProt query expression |
| include_isoforms | no | true or false; defaults to false |

Ready-to-edit examples are in assets/build_manifest.tsv. The test manifest
contains both a species Swiss-Prot row and a clade UniProtKB row.

## Run

On the cluster, use the repository's Singularity configuration and a shared
work directory. The pipeline can build all rows in parallel, subject to the
fetch concurrency configured in nextflow.config.

    nextflow run pipelines/protein-db/main.nf \
      -profile slurm \
      -c config/singularity.config \
      -work-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/long-read/protein-db-work \
      --build_manifest /path/to/build_manifest.tsv \
      --outdir /hps/nobackup/flicek/ensembl/genebuild/jackt/long-read/protein-databases \
      -resume

For a quick wiring check:

    nextflow run pipelines/protein-db/main.nf -stub-run -profile test

## Outputs

Each manifest row is published under:

    <outdir>/databases/<build_id>/database.dmnd
    <outdir>/databases/<build_id>/database_proteins.faa
    <outdir>/databases/<build_id>/provenance.tsv
    <outdir>/databases/<build_id>/versions.yml

For long-read DIAMOND QC, pass either the prebuilt database:

    --diamond_database /path/to/databases/felis_catus_swissprot/database.dmnd

or the published FASTA and let that pipeline build its local database:

    --diamond_reference_proteins /path/to/databases/felis_catus_swissprot/database_proteins.faa

The prebuilt .dmnd is preferred when the database has already been built; the
FASTA remains available for auditability and for other protein-QC tools.

## Design notes

The downloader uses the UniProtKB search endpoint with pagination rather than
one large untracked streaming request. Builds are independent Nextflow tasks,
so -resume reuses completed downloads and DIAMOND databases. The pipeline does
not claim a historical UniProt release unless the API reports one; the
reported release and retrieval timestamp are retained in provenance.tsv.
