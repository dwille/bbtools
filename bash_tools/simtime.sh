#!/bin/sh
# simtime.sh
#   Usage: For the given list of simulations, pull the last simulation time and
#           write to file

ROOT_DIR="."
OUTPUT_DIR="output"
FILENAME="part"

SIMS="128-sedi-per
128-sedi-dir
1024-sedi-per"

for s in $SIMS
  do
    ls -t "$ROOT_DIR/$s/$OUTPUT_DIR/part-"* | head -1 >> simtime
  done
