#!/bin/sh

#SBATCH --partition=devel
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=parallelplot
#SBATCH --output=plot.out

# Single run
# Usage: sbatch submit-plot.sh PLOT_UTILITY SIMULATION DIRECTORY TIME
name="./$1 $2 $3"
srun $name

