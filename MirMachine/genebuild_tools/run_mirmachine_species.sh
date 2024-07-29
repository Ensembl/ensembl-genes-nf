#!/usr/bin/env bash

#SBATCH --job-name=mirMachine_$1
#SBATCH --output=logs/output_%j.log
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --partition=production
#SBATCH --mail-user=$USER@ebi.ac.uk
#SBATCH --mail-type=ALL

main=$1
species=$2
accession=$3

nextflow run $main -profile slurm --species $species  --accession $accession