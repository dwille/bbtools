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

## Number of lines in files function (not including header line)
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # Does not include header

## Define structure class
class structtype():
  pass

## Get info
print "Anisotropy Measures plotting utility"
print ""
#root = raw_input("Simulation root: ")
#if not root.endswith('/'):
#  root = root + '/'
#ts = raw_input("Desired start time: ")
#te = raw_input("Desired end time: ")
root = "../sim/"
ts = "500"
te = "501"

## Sort output files and find number within time range
files = sorted_nicely(glob.glob(root + "data-tetrads/tetrad-data-*"))

# Find number of tetrads
nTetrads = file_len(files[0])

# Initialize numpy arrays
det = np.zeros(nTetrads)
var = np.zeros(nTetrads)
trace = np.zeros(nTetrads)
eigVals = np.zeros((nTetrads, 3))

# Enumerate the columns
detCol = (0)
varCol = (1)
traceCol = (2)
gEigValsCols = (3,4,5)

# Loop over all output files, pull data to structures
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('tetrad-data-'):
    ftime = ftime[12:]

  if float(ftime) >= float(ts) and float(ftime) <= float(te):
    ifile = open(fname)
    ifile.readline()

    # If this gets slow, see
    # -- stackoverlfow 8956832 -- python out of memory on large csv file numpy
    # -- "  " 18259393 -- numpy loading csv too slow compared to matlab
    # -- "  " 26482209 -- fastest way to load huge dat into array
    det = np.genfromtxt(fname, skip_header=1, usecols=detCol)
    var = np.genfromtxt(fname, skip_header=1, usecols=varCol)
    trace = np.genfromtxt(fname, skip_header=1, usecols=traceCol)
    eigVals = np.genfromtxt(fname, skip_header=1, usecols=gEigValsCols)

tetradShape = 27*det / (trace*trace*trace)
tetradVar = 1.5*var / (trace*trace)
principalAxes= eigVals / trace[:,None]

# Initialize histogram bins
nBins = float(25)
iEdges = np.linspace(0,1,nBins + 1)
sEdges = np.linspace(-0.25, 2.0, nBins + 1)

weights = np.ones_like(det)/float(len(det))

pdf = plt.figure()

# I1
axI1 = pdf.add_subplot(321)
axI1 = pdf.subplot2grid((2,3), (0,0), colspan=3)
plt.hist(principalAxes[:,0], bins=nBins, weights=weights, normed=False, color="blue", alpha=0.6)
plt.hist(principalAxes[:,1], bins=nBins, weights=weights, normed=False, color="green", alpha=0.6)
plt.hist(principalAxes[:,2], bins=nBins, weights=weights, normed=False, color="red", alpha=0.6)

axI1.set_xlim([0,1])
axI1.set_ylabel('I1')

# I2
axI2 = pdf.add_subplot(323)
plt.hist(principalAxes[:,1], bins=nBins, weights=weights, normed=False)

axI2.set_xlim([0,1])
axI2.set_ylabel('I2')

# I3
axI3 = pdf.add_subplot(325)
plt.hist(principalAxes[:,2], bins=nBins, weights=weights, normed=False)

axI3.set_xlim([0,1])
axI3.set_ylabel('I3')

# shape
axShape = pdf.add_subplot(322)
plt.hist(tetradShape, bins=sEdges, weights=weights, normed=False)

axShape.set_xlim([-0.25, 2])
axShape.set_ylabel('Shape')

# var
axVar = pdf.add_subplot(324)
plt.hist(tetradVar, bins=nBins, weights=weights, normed=False)

axVar.set_xlim([0,1])
axVar.set_ylabel('Var')

# rog
axR = pdf.add_subplot(326)
plt.hist(np.sqrt(trace), bins=nBins, weights=weights, normed=False)

axR.set_ylabel('R')






plt.show()

