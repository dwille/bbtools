#!/bin/sh

#SBATCH --partition=devel,tesla
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=combine
srun $1 $2
