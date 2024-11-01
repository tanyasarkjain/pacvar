#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Account:
#SBATCH --account=fc_xenopus
#
# Partition:
#SBATCH --partition=savio2
#
# Request one node:
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=3:00:00
#
## Command(s) to run (example):

conda activate env_nf

nextflow run main.nf -profile test,singularity --outdir out -resume
