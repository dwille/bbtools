#!/bin/sh
# submit.sh
#
#SBATCH --partition=gpu
#SBATCH --time=01:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=tet-0500-2.0
#SBATCH --output=tetrad.out

srun ./tetrad_init 
