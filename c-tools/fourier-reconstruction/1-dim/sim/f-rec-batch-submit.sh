#!/bin/sh
# phasevel-submit.sh
#
#SBATCH --partition=shared
#SBATCH --time=4:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=partvel
#SBATCH --output=partout

# Multiple runs

srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho2.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho2.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho2.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho2.0

srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho3.3
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho3.3
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho3.3
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho3.3

srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho4.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho4.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho4.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho4.0

srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/500/rho5.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1000/rho5.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/1500/rho5.0
srun ./f-reconstruct /home-1/dwillen3@jhu.edu/scratch/triply_per/2000/rho5.0
