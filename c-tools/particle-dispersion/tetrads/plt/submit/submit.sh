#!/bin/sh

#SBATCH --partition=devel,tesla
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=combine

# $1 -- program to run
# $2 -- simulation directory
# $3 -- tstart (doesnt matter)


srun $1 $2 $3
