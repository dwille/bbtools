#!/bin/bash
# pull_data.sh
# Usage: executes matlab pull script on all files

TOOL_PATH="/home-1/dwillen3@jhu.edu/bbtools/automate/matlab"
TOOL_NAME="pull_all"
DATA_DIR="/home-1/dwillen3@jhu.edu/scratch/sims/rz_glass"

cd $TOOL_PATH

matlab -nodisplay -nodesktop -r "try; addpath $TOOL_PATH; $TOOL_NAME; catch; end; quit" >\
"$DATA_DIR/pulloutput" 2> "$DATA_DIR/pullerror"
