#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import AutoMinorLocator
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
print ""
print " ---- Alignment traces plotting utility ---- "
print ""
#root = raw_input("Simulation root: ")
#if not root.endswith('/'):
#  root = root + '/'
root = "../sim/"
datadir = root + "data-tetrads/"
infoFile = datadir + "info.dat"

## Sort output files and find nTetrads, nFiles
files = sorted_nicely(glob.glob(datadir + "raw-data-*"))
nTetrads = np.genfromtxt(infoFile, skip_header=1, usecols=0)
nFiles = np.genfromtxt(infoFile, skip_header=1, usecols=1)
time = np.zeros(nFiles)

print "      Found " + str(nFiles) + " files."
print "      Found " + str(nTetrads) + " tetrads."
print ""

## Create array to store alignment of each tetrad at each time
g1s1 = np.zeros((nTetrads, nFiles))

## Pull initial s1 eigenvectors
s0 = np.zeros((nTetrads, 3))
s0[:,0] = np.genfromtxt(files[0], skip_header=1, usecols=18)
s0[:,1] = np.genfromtxt(files[0], skip_header=1, usecols=19)
s0[:,2] = np.genfromtxt(files[0], skip_header=1, usecols=20)


## Loop over files, pull g1s1 alignment at each time for each tetrad
gi = np.zeros((nTetrads, 3))
for ff, fname in enumerate(files):
  gi[:,0] = np.genfromtxt(fname, skip_header=1, usecols=6)
  gi[:,1] = np.genfromtxt(fname, skip_header=1, usecols=7)
  gi[:,2] = np.genfromtxt(fname, skip_header=1, usecols=8)
  
  # Scalar product
  for i in np.arange(3):
    g1s1[:, ff] += gi[:,i]*s0[:,i]

  g1s1[:,ff] = np.abs(g1s1[:,ff])
  # Pull time from filename
  ftime = fname.split("/")[-1]
  time[ff] = float(ftime[9:])

# Find mean
g1s1mean = np.mean(g1s1, 0)

Fig = plt.figure(figsize=(12,8))

ax = Fig.add_subplot(111)
ax.plot(time, g1s1mean, '.', linewidth=3)
for i in np.arange(nTetrads):
  ax.plot(time, g1s1[i,:])

ax.set_xlabel('Time')
ax.set_ylabel('Alignment')

plt.show()
