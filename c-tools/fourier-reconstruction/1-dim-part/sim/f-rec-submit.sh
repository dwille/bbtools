#!/bin/sh
# phasevel-submit.sh
#
#SBATCH --partition=shared
#SBATCH --time=2:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=f-rec-05-rho3.3
#SBATCH --output=f-rec.out
#SBATCH --open-mode=append


srun ./f-reconstruct
