#!/bin/sh
#
#SBATCH --partition=tesla,devel
#SBATCH --gres=gpu:1
#SBATCH --job-name=/home/dwille/scratch/triply_per/2000/rho3.3
#SBATCH --output=tetrad.out
#SBATCH --dependency=after:114390
#SBATCH --open-mode=append
srun ./mod_and_submit.sh $1 $2 $3
