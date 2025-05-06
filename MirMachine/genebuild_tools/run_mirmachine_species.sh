#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --partition=production
#SBATCH --mail-type=FAIL

main=$1
species=$2
accession=$3

nextflow run $main -c nextflow.config -profile slurm --species "$species"  --accession $accession