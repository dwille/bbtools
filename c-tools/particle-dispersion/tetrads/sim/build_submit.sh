#!/bin/bash
# build_submit.sh
# Inputs
#   $1 -- hostname to run on
#   $2 -- job name ($base directory name)
#   $3 -- dependency

echo "#!/bin/sh" > submit.sh
echo "#" >> submit.sh

if [ $1 == "guinevere" ]; then
  echo "#SBATCH --partition=devel" >> submit.sh;
elif [ $1 == "arthur" ]; then
  echo "#SBATCH --partition=tesla" >> submit.sh;
elif [ $1 == "either" ]; then
  echo "#SBATCH --partition=tesla,devel" >> submit.sh;
fi
echo "#SBATCH --gres=gpu:1" >> submit.sh
echo "#SBATCH --job-name=$2" >> submit.sh
echo "#SBATCH --output=tetrad.out" >> submit.sh

chmod +x submit.sh
