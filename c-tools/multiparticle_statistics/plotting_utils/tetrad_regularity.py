#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors, gridspec
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
print " ---- Anisotropy Measures plotting utility ---- "
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

print "      Found " + str(nFiles) + " files in time range"
print ""

# Find number of tetrads
nTetrads = file_len(files[0])

# Initialize numpy arrays
data = [ structtype() for i in range(nFiles) ]

for i in range(nFiles):
  data[i].I1 = np.zeros(nTetrads)
  data[i].I2 = np.zeros(nTetrads)
  data[i].I3 = np.zeros(nTetrads)
  data[i].shape = np.zeros(nTetrads)
  data[i].var = np.zeros(nTetrads)
  data[i].R2 = np.zeros(nTetrads)

time = np.zeros(nFiles)

# Enumerate the columns
R2Col = (0)
varCol = (1)
shapeCol = (2)
maxEig = (3)
medEig = (4)
minEig = (5)

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
    data[i].R2 = np.genfromtxt(fname, skip_header=1, usecols=R2Col)
    data[i].var = np.genfromtxt(fname, skip_header=1, usecols=varCol)
    data[i].shape = np.genfromtxt(fname, skip_header=1, usecols=shapeCol)

    data[i].I1 = np.genfromtxt(fname, skip_header=1, usecols=maxEig)
    data[i].I2 = np.genfromtxt(fname, skip_header=1, usecols=medEig)
    data[i].I3 = np.genfromtxt(fname, skip_header=1, usecols=minEig)

    i += 1

time = time - time[0]

# Initialize histogram bins
nBins = float(35)
weights = np.ones(nTetrads)/nTetrads

initial = plt.figure(figsize=(12,8))
initial.suptitle('Tetrad Regularity at Initialization', fontsize=16)
gs = gridspec.GridSpec(2,3)

# Initial Timestep
axEig = initial.add_subplot(gs[0,:])

n0, bins0, tmp0 = plt.hist(data[0].I1, bins=nBins, weights=weights, 
  normed=False, color="blue", alpha=0.6)
n1, bins1, tmp1 = plt.hist(data[0].I2, bins=nBins, weights=weights, 
  normed=False, color="green", alpha=0.6)
n2, bins2, tmp2 = plt.hist(data[0].I3, bins=nBins, weights=weights, 
  normed=False, color="red", alpha=0.6)

axEig.set_xlim([0,1])
axEig.set_ylim([0,0.12])
axEig.set_xlabel('I_k')
axEig.set_ylabel('Probability')
axEig.legend(['I_1', 'I_2', 'I_3'])
axEig.set_xticks(np.linspace(0, 1, 11))
axEig.xaxis.set_minor_locator(AutoMinorLocator())
axEig.yaxis.set_minor_locator(AutoMinorLocator())
axEig.tick_params(which='major', length=6)
axEig.tick_params(which='minor', length=3)

x0 = bins0[np.argmax(n0)]
y0 = np.max(n0)
axEig.plot(x0, y0, 'o')
axEig.text(x0, y0 + 0.005, '(%.2f' % x0 + ', %.2f' % y0 + ')')

x1 = bins1[np.argmax(n1)]
y1 = np.max(n1)
axEig.plot(x1, y1, 'o')
axEig.text(x1, y1 + 0.005, '(%.2f' % x1 + ', %.2f' % y1 + ')')

x2 = bins2[np.argmax(n2)]
y2 = np.max(n2)
axEig.plot(x2, y2, 'o')
axEig.text(x2, y2 + 0.005, '(%.2f' % x2 + ', %.2f' % y2 + ')')

# shape
axShape = plt.subplot(gs[1,0])
plt.hist(data[0].shape, weights=weights, normed=False, color='black', alpha=0.6)

axShape.set_xlabel('Shape')
axShape.set_ylabel('Probability')
axShape.locator_params(nbins=6)
axShape.xaxis.set_minor_locator(AutoMinorLocator())
axShape.yaxis.set_minor_locator(AutoMinorLocator())
axShape.tick_params(which='major', length=6)
axShape.tick_params(which='minor', length=3)

# var
axVar = plt.subplot(gs[1,1])
plt.hist(data[0].var, weights=weights, normed=False, color='black', alpha=0.6)

axVar.set_xlabel('Variance')
axVar.locator_params(nbins=6)
axVar.xaxis.set_minor_locator(AutoMinorLocator())
axVar.yaxis.set_minor_locator(AutoMinorLocator())
axVar.tick_params(which='major', length=6)
axVar.tick_params(which='minor', length=3)

# rog
axR = plt.subplot(gs[1,2])
plt.hist(np.sqrt(data[0].R2), weights=weights, normed=False, 
  color='black', alpha=0.6)

axR.set_xlabel('R_g')
axR.locator_params(nbins=6)
axR.xaxis.set_minor_locator(AutoMinorLocator())
axR.yaxis.set_minor_locator(AutoMinorLocator())
axR.tick_params(which='major', length=6)
axR.tick_params(which='minor', length=3)


# Time evolution of principal axes
timePrAx = plt.figure(figsize=(12,8))
timePrAx.suptitle('Time Evolution of Principal Axes', fontsize=16)

nPlot = 6
color = np.linspace(1, 0.2, num=nPlot)
tstepplot = np.linspace(0,nFiles-1,nPlot)
tstepplot = tstepplot.astype(int)
legText = ['']*np.size(tstepplot)

i1_ax = timePrAx.add_subplot(311)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].I1, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  i1_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)
  legText[i] = str(time[nn])

i1_ax.set_xlim([0,1])
i1_ax.set_ylabel('P(I_1)')

i2_ax = timePrAx.add_subplot(312)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].I2, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  i2_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

i2_ax.set_xlim([0,1])
i2_ax.set_ylabel('P(I_2)')
i2_ax.legend(legText)

i3_ax = timePrAx.add_subplot(313)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].I3, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  i3_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

i3_ax.set_xlim([0,1])
i3_ax.set_ylabel('P(I_3)')
i3_ax.set_xlabel('I_k')

# Time evolution of shape
timeShape = plt.figure(figsize=(12,8))
timeShape.suptitle('Time Evolution of Shape Measures', fontsize=16)

shape_ax = timeShape.add_subplot(313)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].shape, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  shape_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)
  legText[i] = str(time[nn])

shape_ax.set_xlim([-0.25,2])
shape_ax.set_xlabel('shape')
shape_ax.set_ylabel('P(shape)')
shape_ax.legend(legText)

var_ax = timeShape.add_subplot(312)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].var, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  var_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

var_ax.set_xlim([0,1])
var_ax.set_xlabel('var')
var_ax.set_ylabel('P(var)')

rg_ax = timeShape.add_subplot(311)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(np.sqrt(data[nn].R2), weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  rg_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

rg_ax.set_xlabel('R_g')
rg_ax.set_ylabel('P(R_g)')

plt.show()
