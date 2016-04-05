#!/bin/bash
# tetrad-submit-multiple.sh

#SBATCH --partition=devel
#SBATCH --time=0:05:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=05-2.0-mult
#SBATCH --output=tetrad.out

# Some variables
#DEPENDENCY=""
JOB_FILE="./tetrad_init"

# Set up start times for multiple runs
RUNSTART=({500..510..2})
echo "RUNSTART = ${RUNSTART[*]}"
echo ""

IND=0
for RUNS in "${RUNSTART[@]}"
do
  # Fill in the correct starting time
  sed -i "1s/.*/tStart\ ${RUNSTART[$IND]}/" tetrad.config
  ((IND++))

  # Set up submit command
  JOB_CMD="srun"

  # Add dependency if not first loop
  #if [ -n "$DEPENDENCY" ]
  #then
  #  JOB_CMD="$JOB_CMD --dependency afterok:$DEPENDENCY" 
  #fi

  # Append tetrad_init to JOB_CMD
  JOB_CMD="$JOB_CMD $JOB_FILE"
  echo "Run = $RUNS: Running command: $JOB_CMD"
  OUT=`$JOB_CMD`
  echo "Result: $OUT"
  #DEPENDENCY=`echo $OUT | awk '{print $4}'`

done
