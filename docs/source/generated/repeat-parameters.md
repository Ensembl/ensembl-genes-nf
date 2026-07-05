# Ensembl repeat pipeline parameters

Generated from `pipelines/repeat/nextflow_schema.json`.

## Input/output options

Define where the pipeline should find input data, write results, and cache downloaded files.

| Parameter | Type | Default | Required | Description |
| --- | --- | --- | --- | --- |
| `csvFile` | string |  | yes | Path to the input CSV metadata file. |
| `outdir` | string | ./results | yes | Output directory. |
| `tracedir` | string | ./results/pipeline_info | no | Directory for execution reports and trace files. |
| `cacheDir` | string | ./results/cache | no | Directory used by cached fetch steps. |
| `files_latency` | integer | 60 | no | Sleep time in seconds after file-publishing steps. |
| `cleanup` | boolean | true | no | Remove temporary files when supported. |
| `cleanCache` | boolean | false | no | Remove cached downloaded files after the run. |

## Repeat analysis options

Enable repeat annotation steps and configure repeat tools.

| Parameter | Type | Default | Required | Description |
| --- | --- | --- | --- | --- |
| `generate_lib` | boolean | false | no | Check for an existing RepeatModeler library and generate one when it is missing. |
| `run_repeatmasker` | boolean | false | no | Run RepeatMasker using the RepeatModeler library channel. |
| `run_red` | boolean | false | no | Run RED on the fetched genome assembly. |
| `run_dust` | boolean | false | no | Run DustMasker on the fetched genome assembly. |
| `run_trf` | boolean | false | no | Run TRF on the fetched genome assembly. |
| `ncbiBaseUrl` | string | https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/ | no | NCBI datasets API base URL used for genome downloads. |
| `repeats_ftp_base` | string | http://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/ | no | Base URL for pre-computed RepeatModeler libraries. |
| `repeatmasker_path` | string |  | no | Path to the RepeatMasker executable. |
| `builddatabase_path` | string |  | no | Path to the BuildDatabase executable. |
| `engine_repeatmasker` | string | rmblast | no | RepeatMasker search engine. |
| `engine_repeatmodeler` | string | ncbi | no | RepeatModeler search engine. |
| `red_path` | string |  | no | Path to the RED executable. |
| `dust_path` | string |  | no | Path to the DustMasker executable. |
| `trf_path` | string |  | no | Path to the TRF executable. |
| `trf_match_score` | integer | 2 | no | TRF match score. |
| `trf_mismatch_score` | integer | 5 | no | TRF mismatch penalty. |
| `trf_delta` | integer | 7 | no | TRF indel penalty. |
| `trf_pm` | integer | 80 | no | TRF match probability. |
| `trf_pi` | integer | 10 | no | TRF indel probability. |
| `trf_minscore` | integer | 40 | no | TRF minimum alignment score to report. |
| `trf_maxperiod` | integer | 500 | no | TRF maximum period size. |
