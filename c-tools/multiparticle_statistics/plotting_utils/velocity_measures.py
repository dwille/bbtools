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

## Alignment of two vectors -- <(a.b)^2>
def alignVects(vec1, vec2):
  tmp = np.multiply(vec1, vec2)
  tmp = np.sum(tmp, axis=1)
  tmp = np.square(tmp)
  tmp = np.mean(tmp)
  return tmp

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
te = "2000"

## Sort output files and find number within time range
files = sorted_nicely(glob.glob(root + "data-tetrads/tetrad-data-*"))
inRange = np.zeros(len(files))
n = 0
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('tetrad-data-'):
    time = ftime[12:]

  if float(time) >= float(ts) and float(time) <= float(te):
    inRange[n] = 1
    n += 1

nFiles = n

print "Found " + str(nFiles) + " files in time range"
print ""

# Find number of tetrads
nTetrads = file_len(files[0])

# Initialize numpy arrays
gEigVals = [ structtype() for i in range(nFiles) ]
gEigVecs = [ structtype() for i in range(nFiles) ]
sEigVals = [ structtype() for i in range(nFiles) ]
sEigVecs = [ structtype() for i in range(nFiles) ]
vorticity = [ structtype() for i in range(nFiles) ]
alignment = [ structtype() for i in range(1) ]

for i in range(nFiles):
  gEigVals[i].vals = np.zeros((nTetrads, 3))
  gEigVecs[i].vals = np.zeros((nTetrads, 9))
  sEigVals[i].vals = np.zeros((nTetrads, 3))
  sEigVecs[i].vals = np.zeros((nTetrads, 9))
  vorticity[i].vals = np.zeros((nTetrads, 3))
  vorticity[i].norm = np.zeros((nTetrads, 1))

alignment[0].g1_s1 = np.zeros(nFiles)
alignment[0].g1_s2 = np.zeros(nFiles)
alignment[0].g1_s3 = np.zeros(nFiles)
alignment[0].g2_s1 = np.zeros(nFiles)
alignment[0].g2_s2 = np.zeros(nFiles)
alignment[0].g2_s3 = np.zeros(nFiles)
alignment[0].g3_s1 = np.zeros(nFiles)
alignment[0].g3_s2 = np.zeros(nFiles)
alignment[0].g3_s3 = np.zeros(nFiles)

time = np.zeros(nFiles)

# Enumerate the columns
gEigValsCols = (3,4,5)
gEigVecsCols = (6,7,8,9,10,11,12,13,14)
sEigValsCols = (15,16,17)
sEigVecsCols = (18,19,20,21,22,23,24,25,26)
vorticityCols = (27,28,29)

maxAx = (0,1,2)
medAx = (3,4,5)
minAx = (6,7,8)

# Initialize histogram bins
nBins = float(100)
g1s1 = np.zeros([nBins, nFiles])
edges = np.linspace(0,0.3,nBins + 1)
centers = edges[0:nBins] + 0.5*0.3/nBins

# Loop over all output files, pull data to structures
i = 0;
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('tetrad-data-'):
    ftime = ftime[12:]

  if float(ftime) >= float(ts) and float(ftime) <= float(te):
    time[i] = float(ftime)
    ifile = open(fname)
    ifile.readline()

    # If this gets slow, see
    # -- stackoverlfow 8956832 -- python out of memory on large csv file numpy
    # -- "  " 18259393 -- numpy loading csv too slow compared to matlab
    # -- "  " 26482209 -- fastest way to load huge dat into array

    gEigVals[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      gEigValsCols)
    gEigVecs[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      gEigVecsCols)
    sEigVals[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      sEigValsCols)
    sEigVecs[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      sEigVecsCols)
    vorticity[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      vorticityCols)

    vorticity[i].norm = np.sqrt(np.square(vorticity[i].vals[:,0]) 
                              + np.square(vorticity[i].vals[:,1])
                              + np.square(vorticity[i].vals[:,2]))
    vorticity[i].vals = np.divide(vorticity[i].vals, vorticity[i].norm[:,None])

    # Alignment of principal axes of tetrads with prncpl axes of strain
    alignment[0].g1_s1[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                       sEigVecs[0].vals[:,maxAx])
    alignment[0].g1_s2[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                       sEigVecs[0].vals[:,medAx])
    alignment[0].g1_s3[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                       sEigVecs[0].vals[:,minAx])

    alignment[0].g2_s1[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                       sEigVecs[0].vals[:,maxAx])
    alignment[0].g2_s2[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                       sEigVecs[0].vals[:,medAx])
    alignment[0].g2_s3[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                       sEigVecs[0].vals[:,minAx])

    alignment[0].g3_s1[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                       sEigVecs[0].vals[:,maxAx])
    alignment[0].g3_s2[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                       sEigVecs[0].vals[:,medAx])
    alignment[0].g3_s3[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                       sEigVecs[0].vals[:,minAx])


    # Histogram that shit
    #g1s1[:,i], tmp = np.histogram(alignment[0].g1_s1[:,i], bins=edges,
    #  range=[0,0.3], density=True)

    i += 1


pdf = plt.figure()

g1s_Ax = pdf.add_subplot(311)
g1s_Ax.plot(time, alignment[0].g1_s1[:], 'ko-')
g1s_Ax.plot(time, alignment[0].g1_s2[:], 'bo-')
g1s_Ax.plot(time, alignment[0].g1_s3[:], 'ro-')

g1s_Ax.set_ylabel('g1.s')

g2s_Ax = pdf.add_subplot(312)
g2s_Ax.plot(time, alignment[0].g2_s1[:], 'ko-')
g2s_Ax.plot(time, alignment[0].g2_s2[:], 'bo-')
g2s_Ax.plot(time, alignment[0].g2_s3[:], 'ro-')

g2s_Ax.set_ylabel('g2.s')

g3s_Ax = pdf.add_subplot(313)
g3s_Ax.plot(time, alignment[0].g3_s1[:], 'ko-')
g3s_Ax.plot(time, alignment[0].g3_s2[:], 'bo-')
g3s_Ax.plot(time, alignment[0].g3_s3[:], 'ro-')
g3s_Ax.set_ylabel('g2.s [ms]')

g3s_Ax.set_xlabel('Time [ms]')

plt.show()

