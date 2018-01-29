#!/bin/bash
# loop_initial_times.sh
# Inputs:
#   HOST -- host to run on
#   BASE_DIR -- top-level simulation directory
#   T_STAT -- statistically stationary time
#   TS_INC -- time increment for starting time
#   T_DUR -- duration of run
#   T_FINAL -- simulation final time

# Commandline-args -- check and assign
if [ "$#" -ne 6 ]; then
  echo "Illegal number of parameters"
  echo "Usage: ./tetrad-submit-multiple.sh BASE_DIR T_STAT T_DUR TS_INC"
  echo "   HOST -- guinevere or arthur"
  echo "   BASE_DIR -- top-level simulation directory"
  echo "   T_STAT -- statistically stationary time"
  echo "   TS_INC -- time increment"
  echo "   T_DUR -- duration of a run"
  echo "   T_FINAL -- last simulation time"
  echo ""
  exit 0
else
  HOST=$1
  BASE_DIR=$2
  T_STAT=$3
  TS_INC=$4
  T_DUR=$5
  T_FINAL=$6
fi

# Config file
TETRAD_CONFIG_FILE="$BASE_DIR/analysis/tetrads/tetrad.config"
# Check if config file exists
if [ ! -f $TETRAD_CONFIG_FILE ]; then
  echo "File $TETRAD_CONFIG_FILE not found!"
  exit 0
fi

# Starting run time
RUN_START=$(eval echo {$T_STAT..$(expr $T_FINAL - $T_DUR)..$TS_INC})

# Loop over all starting times, change config file, and submit job
sed -i "11s/.*/Multiple\ Runs\ 1/" $TETRAD_CONFIG_FILE
IND=0
for RUNS in ${RUN_START}; do
  TS=$RUNS
  TE=$(expr $RUNS + $T_DUR)

  echo "Running $BASE_DIR for $TS to $TE"

  # Build submission file
  ./build_submit.sh $HOST $BASE_DIR
  if [ ! -z $JOBID ]; then
    echo "#SBATCH --dependency=after:$JOBID" >> submit.sh
    echo "#SBATCH --open-mode=append" >> submit.sh
  fi
  echo "srun ./mod_and_submit.sh \$1 \$2 \$3" >> submit.sh

  # Submit file
  SCRIPT_SUBMIT="sbatch submit.sh $BASE_DIR $TS $TE"
  JOBID=$($SCRIPT_SUBMIT)
  JOBID=`echo $JOBID | awk '{print $4}'`

  echo ""
  # Increment to next time series
  ((IND++))
done

# Return config file to initial
# sed -i "1s/.*/tStart\ $T_STAT/" $TETRAD_CONFIG_FILE
# sed -i "2s/.*/tEnd\ 25000/" $TETRAD_CONFIG_FILE
# sed -i "11s/.*/Multiple\ Runs\ 0/" $TETRAD_CONFIG_FILE

# # Print config
# echo "Base simulation directory = $BASE_DIR"
# echo "Tetrad config file = $TETRAD_CONFIG_FILE"
# echo ""
# echo "Initial Starting time = $T_STAT"
# echo "Time increment = $TS_INC"
# echo "Starting times = ${RUN_START[*]}"
# echo "Duration = $T_DUR"
# echo ""

