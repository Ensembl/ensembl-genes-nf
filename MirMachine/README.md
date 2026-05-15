# ensembl-miRNA-nf
Automated MirMachine Based miRNA Detection workflow developed for genebuild team at Ensembl

## Deployment
The only requirements for running this workflow from a software perspective are [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Singularity/Apptainer](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages). 

Original development and testing was on EMBL-EBIs Codon HPC with further testing on a standalone Ubuntu 20.04 machine. Tested with Nextflow version 24.10.4 and Apptainer 1.3.4

### General

Default inputs to this Nextflow pipeline are specified in the params.config file. However, each can be overwritten when calling Nextflow:
Eg.
```
nextflow run main.nf \
  --input <path to tsv> \ # Tab separated file with columns titled 'Scientific Name' and 'Accession'. These must match what is in Ensembl FTP
  --outdir <path to output directory>
  --fasta_dir <path to where fastas are stored> # Used to speed up processing when fastas are on local file system
```

## MIR-12 Lepidoptera analysis

This branch is for the MirMachine 2 manuscript rerun requested in May 2026:
MIR-12 only, across the Lepidoptera genomes supplied by Vanessa Molin Paynter.
It is intentionally analysis-specific. General MirMachine workflow improvements
should stay on `feature/mirmachine`.

The run uses:

- input file: `inputs/lepidoptera_mir12_input.tsv`
- family: `Mir-12`
- model: `combined`
- e-value: `5`
- long hairpin covariance models: enabled with `--long`
- MirMachine source: `sinanugur/MirMachine` development commit
  `14e33308cbeaed495c24254f60bd8c52232d7b2e`

The container image is configured in `params.config` as:

```text
docker://ensemblorg/mirmachine:development-14e3330
```

The Dockerfile for this image belongs in the
`github.com/ensembl/ensembl-genes-containers` repository under
`Containers/mirmachine/`. If the image is published under a different name, pass
the replacement at runtime with `--mirmachine_container`.

Run on Slurm:

```bash
./run_mir12_lepidoptera.sh
```

Override output and FASTA locations:

```bash
MIR12_OUTDIR=/path/to/results MIR12_FASTA_DIR=/path/to/fastas ./run_mir12_lepidoptera.sh
```

This analysis branch was prepared with assistance from Codex. Runtime validation
and result interpretation remain the responsibility of the project authors.
