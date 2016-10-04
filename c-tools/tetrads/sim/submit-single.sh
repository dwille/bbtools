#!/bin/sh
#
# $1 -- dir to run
echo $1
#SBATCH --partition=devel,tesla
#SBATCH --gres=gpu:1
#SBATCH --job-name=tetrads
#SBATCH --output=tetrad.out
#SBATCH --open-mode=append
srun ./tetrads $1
