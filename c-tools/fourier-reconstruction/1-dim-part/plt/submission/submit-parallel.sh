#!/bin/sh
# Submit multiple plotting jobs to slurm, in parallel
# Input args:
# 1:  Plot utility program name

PlotUtility=$1

sbatch submit-single-plot.sh $1 500/rho2.0 400
sbatch submit-single-plot.sh $1 500/rho3.3 400
sbatch submit-single-plot.sh $1 500/rho4.0 400
sbatch submit-single-plot.sh $1 500/rho5.0 400

sbatch submit-single-plot.sh $1 1000/rho2.0 650
sbatch submit-single-plot.sh $1 1000/rho3.3 450
sbatch submit-single-plot.sh $1 1000/rho4.0 400
sbatch submit-single-plot.sh $1 1000/rho5.0 670
                                               
sbatch submit-single-plot.sh $1 1500/rho2.0 800
sbatch submit-single-plot.sh $1 1500/rho3.3 400
sbatch submit-single-plot.sh $1 1500/rho4.0 450
                                
sbatch submit-single-plot.sh $1 2000/rho2.0 600
sbatch submit-single-plot.sh $1 2000/rho3.3 550

echo "Finished" >> plot.out
