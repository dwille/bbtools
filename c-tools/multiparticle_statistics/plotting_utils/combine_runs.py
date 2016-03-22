#!/usr/bin/env python2
# Combine various runs of a tetrad analysis into master dat files

import sys, os, re
import glob
import numpy as np
import matplotlib.pyplot as plt

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line

## Nicely sorted function
def sorted_nicely( l ):
  """ Sorts the given iterable in a natural way

  Required arguments:
  l -- The iterable to be sroted.

  courtesy stackoverflow/2669059
  """
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
  return sorted(l, key = alphanum_key)

print "--- Combine-Runs Utility ---"
print ""

# Set root dir
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "  Sim root directory set to: " + root

# Get simdir and create datadir
simdir = raw_input("  Simulation directory: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'
datadir = root + simdir + "data-tetrads/"

# Get list of directories in the datadir -- corresponds to different runs
runs = sorted_nicely(os.listdir(datadir))

# Initalize data structs and arrays
nTSteps = np.zeros(len(runs))
data = [ structtype() for i in range(len(runs)) ]

# Loop over all runs, find nTSteps
for rr, run in enumerate(runs):
  rundir = datadir + run + '/'
  statFile = rundir + 'stat.dat'
  nodeFile = rundir + 'regularNodes'

  nTSteps[rr] = file_len(statFile)

# find maxTSteps
maxTSteps = np.max(nTSteps)
print maxTSteps

# init numerator, denominator arrays
numR2 = np.zeros(maxTSteps)
numVar = np.zeros(maxTSteps)
numShape = np.zeros(maxTSteps)
denom = np.zeros(maxTSteps)

for rr, run in enumerate(runs):
  rundir = datadir + run + '/'
  statFile = rundir + 'stat.dat'
  alignFile = rundir + 'align.dat'
  nodeFile = rundir + 'regularNodes'

  # init data to maxTStep size
  data[rr].nTetrads = np.zeros(maxTSteps)
  data[rr].time = np.zeros(maxTSteps)
  data[rr].meanR2 = np.zeros(maxTSteps)
  data[rr].meanVar = np.zeros(maxTSteps)
  data[rr].meanShape = np.zeros(maxTSteps)

  # only fill up to nTSteps[rr] -- this pads the end of the array with 0s
  time = np.genfromtxt(statFile, skip_header=1, usecols=0)
  data[rr].time[0:nTSteps[rr]] = time - time[0]
  data[rr].meanR2[0:nTSteps[rr]] = np.genfromtxt(statFile, skip_header=1, 
    usecols=1)
  data[rr].meanVar[0:nTSteps[rr]] = np.genfromtxt(statFile, skip_header=1, 
    usecols=2)
  data[rr].meanShape[0:nTSteps[rr]] = np.genfromtxt(statFile, skip_header=1, 
    usecols=3)
  data[rr].nTetrads[0:nTSteps[rr]] = file_len(nodeFile)*np.ones(nTSteps[rr])

  numR2 += data[rr].meanR2*data[rr].nTetrads
  numVar += data[rr].meanR2*data[rr].nTetrads
  numShape += data[rr].meanR2*data[rr].nTetrads
  denom += data[rr].nTetrads

# Find overall mean
print 'Total tetrads tracked: ' + str(denom[0])

meanR2 = numR2 / denom
meanVar = numVar / denom
meanShape = numShape / denom

for rr in np.arange(len(runs)):
  plt.plot(data[rr].time)

plt.show()


