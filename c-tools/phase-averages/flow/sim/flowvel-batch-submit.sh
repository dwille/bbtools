#!/bin/sh
# phasevel-submit.sh
#
#SBATCH --partition=debug
#SBATCH --time=5:0:0
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --job-name=batch-flowvel
#SBATCH --output=flowphaseout

# Multiple runs

srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho2.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho2.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho2.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho2.0

srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho3.3
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho3.3
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho3.3
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho3.3

srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho4.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho4.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho4.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho4.0

srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho5.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho5.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho5.0
srun ./flowVel /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho5.0
