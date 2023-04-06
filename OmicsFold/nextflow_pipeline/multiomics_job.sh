#!/bin/bash

#script to submit nextflow multiomics job to slurm job scheduler

#SBATCH -J multiomics_nextflow
#SBATCH -o multiomics_%j_output.txt
#SBATCH -e multiomics_%j_error.txt
#SBATCH -N 1
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128GB

nextflow run multiomics_nextflow.nf.groovy --data $1 --data_labels $2

#example usage on command line
#sbatch multiomics_job.sh data.rds data_labels.rds
