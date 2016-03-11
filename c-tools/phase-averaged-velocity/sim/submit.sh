#!/bin/sh
# submit.sh
#
#SBATCH --partition=devel
#SBATCH --gres=gpu:1
#SBATCH --job-name=test-tetrad
#SBATCH --output=out

#srun valgrind --leak-check=full\
#              --suppressions=sppressions\
#              ./tetrad_init ./output/part-0.0.cgns 2.5 .5

#srun valgrind --gen-suppressions=all\
#              --leak-check=full\
#              --suppressions=sppressions\
#              --show-reachable=yes\
#              ./tetrad_init ./output/part-0.0.cgns 2.5 .5

#srun cuda-memcheck --leak-check full\
#                   ./tetrad_init ./output/part-0.0.cgns 2.5 .5

srun ./main
