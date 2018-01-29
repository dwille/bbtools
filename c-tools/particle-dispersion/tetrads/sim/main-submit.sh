#!/bin/sh
# Submit multiple plotting jobs to slurm, in parallel
# This is the first script to run!
# Input args:
# 1:  Plot utility program name

# loop_initial_times.sh
#       HOST
#       BASE_DIR
#       T_STAT
#       TS_INC
#       T_DUR
#       T_FINAL

./loop_initial_times.sh either /home/dwille/scratch/triply_per/500/rho2.0 400 300 6000 25000
./loop_initial_times.sh either /home/dwille/scratch/triply_per/500/rho3.3 400 154 4180 21336
./loop_initial_times.sh either /home/dwille/scratch/triply_per/500/rho4.0 400 160 1952 21532
./loop_initial_times.sh either /home/dwille/scratch/triply_per/500/rho5.0 400 116 1892 18498

./loop_initial_times.sh either /home/dwille/scratch/triply_per/1000/rho2.0 650 264 6662 25000
./loop_initial_times.sh either /home/dwille/scratch/triply_per/1000/rho3.3 450 174 3952 20410
./loop_initial_times.sh either /home/dwille/scratch/triply_per/1000/rho4.0 400 136 2480 17150
./loop_initial_times.sh either /home/dwille/scratch/triply_per/1000/rho5.0 670 116 2430 16008

./loop_initial_times.sh either /home/dwille/scratch/triply_per/1500/rho2.0 800 362 7726 22216
./loop_initial_times.sh either /home/dwille/scratch/triply_per/1500/rho3.3 400 188 4280 15792
./loop_initial_times.sh either /home/dwille/scratch/triply_per/1500/rho4.0 450 160 3504 14618

./loop_initial_times.sh either /home/dwille/scratch/triply_per/2000/rho2.0 600 564 13298 19370
./loop_initial_times.sh either /home/dwille/scratch/triply_per/2000/rho3.3 550 304 6594 15210
