#!/bin/sh
# phasevel-submit.sh
#
#SBATCH --partition=tesla
#SBATCH --time=6:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=f-rec
#SBATCH --output=rec-log

# Multiple runs

#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/500/rho2.0
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1000/rho2.0
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1500/rho2.0
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/2000/rho2.0
#
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/500/rho3.3
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1000/rho3.3
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1500/rho3.3
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/2000/rho3.3
#
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/500/rho4.0
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1000/rho4.0
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1500/rho4.0
##srun ./f-rec-1D-part /home/dwille/scratch/triply_per/2000/rho4.0
#
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/500/rho5.0
#srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1000/rho5.0
##srun ./f-rec-1D-part /home/dwille/scratch/triply_per/1500/rho5.0
##srun ./f-rec-1D-part /home/dwille/scratch/triply_per/2000/rho5.0

srun ../plt/histogram_vfrac_vel.py 500/rho2.0 0
srun ../plt/histogram_vfrac_vel.py 500/rho3.3 0
srun ../plt/histogram_vfrac_vel.py 500/rho4.0 0
srun ../plt/histogram_vfrac_vel.py 500/rho5.0 0

srun ../plt/histogram_vfrac_vel.py 1000/rho2.0 0
srun ../plt/histogram_vfrac_vel.py 1000/rho3.3 0
srun ../plt/histogram_vfrac_vel.py 1000/rho4.0 0
srun ../plt/histogram_vfrac_vel.py 1000/rho5.0 0

srun ../plt/histogram_vfrac_vel.py 1500/rho2.0 0
srun ../plt/histogram_vfrac_vel.py 1500/rho3.3 0
srun ../plt/histogram_vfrac_vel.py 1500/rho4.0 0

srun ../plt/histogram_vfrac_vel.py 2000/rho2.0 0
srun ../plt/histogram_vfrac_vel.py 2000/rho3.3 0
