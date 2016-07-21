#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import re

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

## Number of lines in files function (not includinng header line)
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # Does not include header

## safe ln
def safe_ln(x, minval=1e-12):
  return np.log(x.clip(min=minval))

## Get info
print "Anisotropy Measures plotting utility"
print ""
#root = raw_input("Simulation root: ")
#if not root.endswith('/'):
#  root = root + '/'
#ts = raw_input("Desired start time: ")
#te = raw_input("Desired end time: ")
root = "../sim/"
ts = "0"
te = "2000"

## Sort output files and find number within time range
files = sorted_nicely(glob.glob(root + "data-tetrads/raw-data-*"))
inRange = np.zeros(len(files))
n = 0
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('raw-data-'):
    time = ftime[9:]

  if float(time) >= float(ts) and float(time) <= float(te):
    inRange[n] = 1
    n += 1

nFiles = n

print "Found " + str(nFiles) + " files in time range"
print ""

# Find number of tetrads
nTetrads = file_len(files[0])

maxEig = np.zeros((nTetrads, 1))
medEig = np.zeros((nTetrads, 1))
minEig = np.zeros((nTetrads, 1))
avgMaxEig = np.zeros(nFiles)
avgMedEig = np.zeros(nFiles)
avgMinEig = np.zeros(nFiles)
lMultMax = np.zeros(nFiles)
lMultMed = np.zeros(nFiles)
lMultMin = np.zeros(nFiles)
time = np.zeros(nFiles)

# Loop over all output files, pull data to structures
i = 0;
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('raw-data-'):
    ftime = ftime[9:]

  if float(ftime) >= float(ts) and float(ftime) <= float(te):
    time[i] = float(ftime)
    ifile = open(fname)
    ifile.readline()

    # If this gets slow, see
    # -- stackoverlfow 8956832 -- python out of memory on large csv file numpy
    # -- "  " 18259393 -- numpy loading csv too slow compared to matlab
    # -- "  " 26482209 -- fastest way to load huge dat into array

    maxEig = np.genfromtxt(fname, skip_header=1, usecols=3)
    medEig = np.genfromtxt(fname, skip_header=1, usecols=4)
    minEig = np.genfromtxt(fname, skip_header=1, usecols=5)

    avgMaxEig[i] = np.mean(maxEig)
    avgMedEig[i] = np.mean(medEig)
    avgMinEig[i] = np.mean(minEig)

    lMultMax[i] = np.mean(0.5*safe_ln(maxEig)/time[i])
    lMultMed[i] = np.mean(0.5*safe_ln(medEig)/time[i])
    lMultMin[i] = np.mean(0.5*safe_ln(minEig)/time[i])

    i += 1
  #
#

mult = plt.figure()

regAx = mult.add_subplot(311)
regAx.plot(time, avgMaxEig)
regAx.plot(time, avgMedEig)
regAx.plot(time, avgMinEig)
regAx.set_ylabel('avgEig')

lnAx = mult.add_subplot(312)
lnAx.plot(time, safe_ln(avgMaxEig))
lnAx.plot(time, safe_ln(avgMedEig))
lnAx.plot(time, safe_ln(avgMinEig))
lnAx.set_ylabel('ln(avgMaxEig)')

multAx = mult.add_subplot(313)
multAx.plot(time, lMultMax)
multAx.plot(time, lMultMed)
multAx.plot(time, lMultMin)
multAx.set_ylabel('ln(avgMaxEig)/2*time')

plt.show()
