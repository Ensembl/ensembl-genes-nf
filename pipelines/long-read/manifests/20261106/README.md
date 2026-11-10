# Feature/long-read manifests (2026-11-06)

These manifests were generated from the seven-species ENA long-read inventory and the joined BioSamples run metadata. The `feature/long-read` pipeline requires local FASTQ paths; it does not consume the remote `source_uri` directly.

Each species directory contains:

- `primary.tsv`: directly usable after the local `fastq/` files have been staged and MD5-checked. It excludes runs without FASTQ/MD5, amplicon/clone/poolclone strategies, unusual `OTHER`/`ncRNA-Seq` strategies, direct-RNA title matches, raw-subread filename matches, and multi-FASTQ runs.
- `review.tsv`: usable artifacts requiring an explicit decision before processing. This includes raw PacBio subreads, unusual strategies, direct-RNA title matches, and multi-FASTQ runs.
- `quarantine.tsv`: runs excluded from the current workflow because there is no normalized FASTQ+MD5 or the library strategy is amplicon/clone/poolclone.

The `fastq` paths are relative to the repository launch directory, matching the `feature/long-read` pipeline contract. When launching from the feature worktree root, download each `source_uri` into the exact path under `pipelines/long-read/manifests/20261106/<species>/fastq/` and verify `selected_artifact_md5` before launching Nextflow. Multi-FASTQ review rows require concatenation or a pipeline-specific split decision; they are not silently reduced to one file.

The primary manifest contract is `run_id`, `fastq`, `technology`, and `minimap2_preset`; the extra columns preserve provenance and biology. The primary manifests contain no raw PacBio subread rows, so they are the conservative starting point for the feature pipeline.
