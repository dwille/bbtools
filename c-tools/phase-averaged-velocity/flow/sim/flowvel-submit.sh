#!/bin/sh
# phasevel-submit.sh
#
#SBATCH --partition=gpu
#SBATCH --time=0:6:0
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=phasevel
#SBATCH --output=phaseout
#SBATCH --open-mode=append


srun ./phasevel
