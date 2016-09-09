#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import re

## Define structure class
class structtype():
  pass

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

## Get info
print ""
print " ---- Compressibility Plotting Utility ----"
print ""

# Set root dir
#root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
root = "/home/dwille/bbtools/c-tools/multiparticle_statistics/"
print "      Sim root directory set to: " + root + " ..."
print ""

#ts = raw_input("Desired start time: ")
#te = raw_input("Desired end time: ")
ts = "500"
te = "1000"
simdir = "sim/"
datadir = root + simdir + "data-tetrads/"
infoFile = datadir + "info.dat"

## Sort output files and find number within time range
files = sorted_nicely(glob.glob(datadir + "raw-data-*"))
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

print "      Found " + str(nFiles) + " files in time range..."
print ""

# Find number of tetrads
nTetrads = np.genfromtxt(infoFile, skip_header=1, usecols=0)
nTsteps = np.genfromtxt(infoFile, skip_header=1, usecols=1)

# Initialize data structs
RoG = np.zeros((nTetrads, nTsteps))
EVar = np.zeros((nTetrads, nTsteps))
Shape = np.zeros((nTetrads, nTsteps))
S11 = np.zeros((nTetrads, nTsteps))
S22 = np.zeros((nTetrads, nTsteps))
S33 = np.zeros((nTetrads, nTsteps))
time = np.zeros(nTsteps)

# Loop over all output files, pull data to structures
i = 0;
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('raw-data-'):
    ftime = ftime[9:]

  if float(ftime) >= float(ts) and float(ftime) <= float(te):
    time[i] = float(ftime)

    RoG[:,i] = np.genfromtxt(fname, skip_header=1, usecols=0);
    EVar[:,i] = np.genfromtxt(fname, skip_header=1, usecols=1);
    Shape[:,i] = np.genfromtxt(fname, skip_header=1, usecols=2);

    S11[:,i] = np.genfromtxt(fname, skip_header=1, usecols=30);
    S22[:,i] = np.genfromtxt(fname, skip_header=1, usecols=31);
    S33[:,i] = np.genfromtxt(fname, skip_header=1, usecols=32);


    i += 1
  #
#


## Strain vs R
sr_fig = plt.figure()
sr_fig.suptitle('Strain vs R')

s11r_Ax = sr_fig.add_subplot(311)
for i in np.arange(0, nTsteps):
  s11r_Ax.plot(RoG[:,i], S11[:,i], '.')
s11r_Ax.set_ylim([-5,5])
s11r_Ax.set_xlabel('R')
s11r_Ax.set_ylabel('S11')

s22r_Ax = sr_fig.add_subplot(312)
s22r_Ax.plot(RoG[:,:], S22[:,:], '.')
s22r_Ax.set_ylim([-5,5])
s22r_Ax.set_xlabel('R')
s22r_Ax.set_ylabel('S22')

s33r_Ax = sr_fig.add_subplot(313)
s33r_Ax.plot(RoG[:,:], S33[:,:], '.')
s33r_Ax.set_ylim([-5,5])
s33r_Ax.set_xlabel('R')
s33r_Ax.set_ylabel('S33')

## Strain vs Shape
ss_fig = plt.figure()
ss_fig.suptitle('Strain vs Shape')

s11s_Ax = ss_fig.add_subplot(311)
s11s_Ax.plot(Shape[:,:], S11[:,:], '.')
s11s_Ax.set_ylim([-1,1])
s11s_Ax.set_xlabel('Shape')
s11s_Ax.set_ylabel('S11')

s22s_Ax = ss_fig.add_subplot(312)
s22s_Ax.plot(Shape[:,:], S22[:,:], '.')
s22s_Ax.set_ylim([-1,1])
s22s_Ax.set_xlabel('Shape')
s22s_Ax.set_ylabel('S22')

s33s_Ax = ss_fig.add_subplot(313)
s33s_Ax.plot(Shape[:,:], S33[:,:], '.')
s33s_Ax.set_ylim([-1,1])
s33s_Ax.set_xlabel('Shape')
s33s_Ax.set_ylabel('S33')

## Strain vs EVar
sev_fig = plt.figure()
sev_fig.suptitle('Strain vs EVar')

s11ev_Ax = sev_fig.add_subplot(311)
s11ev_Ax.plot(EVar[:,:], S11[:,:], '.')
s11ev_Ax.set_ylim([-1,1])
s11ev_Ax.set_xlabel('EVar')
s11ev_Ax.set_ylabel('S11')

s22ev_Ax = sev_fig.add_subplot(312)
s22ev_Ax.plot(EVar[:,:], S22[:,:], '.')
s22ev_Ax.set_ylim([-1,1])
s22ev_Ax.set_xlabel('EVar')
s22ev_Ax.set_ylabel('S22')

s33ev_Ax = sev_fig.add_subplot(313)
s33ev_Ax.plot(EVar[:,:], S33[:,:], '.')
s33ev_Ax.set_ylim([-1,1])
s33ev_Ax.set_xlabel('EVar')
s33ev_Ax.set_ylabel('S33')

plt.show()

