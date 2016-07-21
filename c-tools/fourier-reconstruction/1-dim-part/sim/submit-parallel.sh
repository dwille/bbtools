#!/bin/sh
# Submit multiple plotting jobs to slurm, in parallel
# Input args:
# 1:  Plot utility program name

sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/500/rho2.0
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/500/rho3.3
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/500/rho4.0
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/500/rho5.0

sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1000/rho2.0
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1000/rho3.3
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1000/rho4.0
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1000/rho5.0

sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1500/rho2.0
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1500/rho3.3
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/1500/rho4.0

sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/2000/rho2.0
sbatch f-rec-submit.sh /home/dwille/scratch/triply_per/2000/rho3.3
