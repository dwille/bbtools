#!/bin/sh
# Submit multiple f-reconstructions at once
#
#SBATCH --partition=devel
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=all
#SBATCH --output=allout


srun ./f-reconstruct $1
