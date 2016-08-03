#!/bin/sh
# simtime.sh
#   Usage: For the given list of simulations, pull the last simulation time and
#           write to file

ROOT_DIR="/home-1/dwillen3@jhu.edu/scratch/sims/rz_glass"
OUTPUT_DIR="output"
FILENAME="part"
OUTPUT_FILE="lastSimTime"

SIMS="fluid_bot/565/dir-020
fluid_bot/565/dir-060
fluid_bot/565/dir-100
fluid_bot/565/per-020
fluid_bot/565/per-060
fluid_bot/565/per-100
fluid_bot/3045/dir-020
fluid_bot/3045/dir-040
fluid_bot/3045/dir-060
fluid_bot/3045/dir-080
fluid_bot/3045/dir-100
fluid_bot/3045/per-020
fluid_bot/3045/per-040
fluid_bot/3045/per-060
fluid_bot/3045/per-080
fluid_bot/3045/per-100
fluid_triply_per/565/565_rho2.0
fluid_triply_per/565/565_rho3.3
fluid_triply_per/565/565_rho4.0
fluid_triply_per/565/565_rho5.0
fluid_triply_per/3045/3045_rho2.0
fluid_triply_per/3045/3045_rho3.3
fluid_triply_per/3045/3045_rho4.0
fluid_triply_per/3045/3045_rho5.0"

if [ -f $ROOT_DIR/$OUTPUT_FILE ];
then
  rm "$ROOT_DIR/$OUTPUT_FILE"
fi

cd $ROOT_DIR

for s in $SIMS
  do
    ls -t "$ROOT_DIR/$s/$OUTPUT_DIR/$FILENAME-"* | head -1 >> "$ROOT_DIR/$OUTPUT_FILE"
  done
