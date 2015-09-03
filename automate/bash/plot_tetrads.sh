#!/bin/bash
# plot_tetrads.sh
# Usage: executes matlab plot script on all files

TOOL_PATH="/home-1/dwillen3@jhu.edu/bbtools/automate/matlab"
TOOL_NAME="plot_tetrads"
DATA_DIR="/home-1/dwillen3@jhu.edu/scratch/sims/rz_glass"

cd $TOOL_PATH

matlab -nodisplay -nodesktop -r "try; addpath $TOOL_PATH; $TOOL_NAME; catch; end; quit" >\
  "$DATA_DIR/plotoutput" 2> "$DATA_DIR/ploterror"
