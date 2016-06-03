#!/bin/sh
# partvel-submit.sh
#
#SBATCH --partition=shared
#SBATCH --time=2:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=partvel
#SBATCH --output=partout

srun ./partVel
