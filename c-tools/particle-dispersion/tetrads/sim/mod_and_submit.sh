#!/bin/bash
# mod_and_submit.sh
# Inputs
#   BASE_DIR -- top-level simulatoin directory
#   TS -- Tetrad sim start time
#   TE -- Tetrad sim end time

# Parse input args
if [ "$#" -ne 3 ]; then
  echo "$#"
  echo "Illegal number of parameters in submit.sh"
  echo "Usage: ./mod_and_submit.sh BASE_DIR TS TE"
  echo ""
  exit 0
else
  BASE_DIR=$1
  TS=$2
  TE=$3
fi

echo "Running $1 from $2 to $3 at $("date")"

# Set input file
TETRAD_CONFIG_FILE="$BASE_DIR/analysis/tetrads/tetrad.config"
if [ ! -f $TETRAD_CONFIG_FILE ]; then
  echo "File $TETRAD_CONFIG_FILE not found!"
  exit 0
fi

# Modify input file
sed -i "1s/.*/tStart\ $TS/" $TETRAD_CONFIG_FILE
sed -i "2s/.*/tEnd\ $TE/" $TETRAD_CONFIG_FILE

# Submit the binary file
./tetrads $BASE_DIR
echo ""
