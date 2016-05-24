#!/bin/sh

#SBATCH --partition=parallel
#SBATCH --time=2:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=plots
#SBATCH --output=plot.out

# Multiple run
# Usage: sbatch submit-plot.sh PLOT_UTILITY
n0500r20="./$1 500/rho2.0 400"
n0500r33="./$1 500/rho3.3 400"
n0500r40="./$1 500/rho4.0 400"
n0500r50="./$1 500/rho5.0 400"

n1000r20="./$1 1000/rho2.0 650"
n1000r33="./$1 1000/rho3.3 450"
n1000r40="./$1 1000/rho4.0 400"
n1000r50="./$1 1000/rho5.0 670"

n1500r20="./$1 1500/rho2.0 800"
n1500r33="./$1 1500/rho3.3 400"
n1500r40="./$1 1500/rho4.0 450"

n2000r20="./$1 2000/rho2.0 600"
n2000r33="./$1 2000/rho3.3 550"

srun $n0500r20
srun $n0500r33
srun $n0500r40
srun $n0500r50

srun $n1000r20
srun $n1000r33
srun $n1000r40
srun $n1000r50

srun $n1500r20
srun $n1500r33
srun $n1500r40

srun $n2000r20
srun $n2000r33

n1500r50="./$1 1500/rho5.0 600"
n2000r40="./$1 2000/rho4.0 600"
n2000r50="./$1 2000/rho5.0 500"
#srun $n1500r50
#srun $n2000r40
#srun $n2000r50
