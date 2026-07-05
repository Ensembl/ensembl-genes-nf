# RiboSeq Pipeline Parameters

Generated from `pipelines/riboseq/nextflow_schema.json`.

## Parameters

Parameter schema for the ribosome profiling pipeline

| Parameter | Type | Default | Required | Description |
| --- | --- | --- | --- | --- |
| `sample_sheet` | string |  | yes | Path to CSV sample sheet with Run and study_accession columns |
| `outdir` | string | ./results | no | Output directory for pipeline results |
| `star_index` | string |  | yes | Path to STAR genome index directory |
| `gtf` | string |  | yes | Path to GTF annotation file |
| `fasta` | string |  | yes | Path to reference genome FASTA file (required for RiboWaltz P-site analysis) |
| `ribometric_annotation` | string |  | no | Path to RiboMetric annotation file (optional, for RiboMetric QC) |
| `chrom_sizes_file` | string |  | yes | Path to chromosome sizes file |
| `adapter_list` | string | ${projectDir}/resources/adapter_list.tsv | no | Path to adapter list TSV file |
| `rrna_index` | string |  | no | Path to bowtie1 rRNA index directory (e.g., /path/to/rRNA_index) or path to any index file. All *.ebwt files in the directory will be staged for rRNA filtering. |
| `fetch` | boolean | true | no | Enable data fetching from SRA |
| `force_fetch` | boolean | false | no | Force reprocessing of existing collapsed reads |
| `collapsed_read_path` | string |  | no | Path to directory containing pre-collapsed reads |
| `use_architecture_detection` | boolean | true | no | Use getRPF architecture detection for adapter finding |
| `mismatches` | integer | 3 | no | Maximum number of mismatches allowed in alignment |
| `alignment_type` | string; one of: EndToEnd, Local | EndToEnd | no | STAR alignment type |
| `allow_introns` | boolean | true | no | Allow introns in STAR alignment |
| `max_multimappers` | integer | 10 | no | Maximum number of multi-mapping locations |
| `min_mapq` | integer | 0 | no | Minimum mapping quality for filtering |
| `trim_front` | integer | 0 | no | Number of bases to trim from 3' end |
| `ribowaltz_exclude_start` | integer | 0 | no | Number of nucleotides to exclude from start codon in RiboWaltz coverage |
| `ribowaltz_exclude_stop` | integer | 0 | no | Number of nucleotides to exclude from stop codon in RiboWaltz coverage |
| `enable_unique_reads_tracking` | boolean | false | no | Create unique-read tracking outputs from collapsed FASTA files |
| `unique_reads_mode` | string; one of: per_run, progressive | per_run | no | Unique-read tracking mode |
| `unique_reads_previous_index` | string |  | no | Path to a previous unique-read index directory when using progressive mode |
