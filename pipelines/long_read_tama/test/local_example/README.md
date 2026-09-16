# Local long-read TAMA example

This example creates a two-read gzipped ONT FASTQ and a three-contig reference,
then runs the pipeline locally in stub mode. Stub mode exercises the complete
workflow without requiring the production TAMA container or a large genome.

Run it from the repository root:

```bash
bash pipelines/long_read_tama/test/local_example/run_local_example.sh
```

The script demonstrates three supported paths:

1. `candidate_manifest.csv` through `--manifest`, including local inventory
   generation and classification reports.
2. The reviewed `approved_manifest.tsv` through `--approved_manifest`, with
   `--shard_mode contig`.
3. The same approved input with explicit `--shard_mode none`.

The candidate CSV is accepted as an input convenience and is normalised to the
pipeline's internal TSV contract. The HTTP server is local-only and allows the
normal acquisition/checksum path to be tested without an archive download.

If `minimap2` and `samtools` are installed, the script also performs a real
alignment and invokes the workload inspector and contig splitter, producing
indexed BAMs under `real-sharding/shards/`. The Nextflow runs use the real
container-backed processes by default. Set `MODE=stub` only for a wiring-only
fallback when Singularity/container execution is unavailable.
