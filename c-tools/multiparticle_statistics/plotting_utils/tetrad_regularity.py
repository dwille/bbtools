#!/usr/bin/env python2

import sys, os
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
print " ---- Tetrad regularity plotting utility ---- "
print ""

# DEVEL
#root = "../sim/"
#datadir = root + "data-tetrads/"

# MARCC
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "      Sim root directory set to: " + root
simdir = raw_input("      Simulation directory: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(datadir):
  print "      " + datadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + simdir + "/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

datadir = root + simdir + "data-tetrads/"
infoFile = datadir + "info.dat"

## Sort output files and find number within time range
files = sorted_nicely(glob.glob(datadir + "raw-data-*"))
nTetrads = np.genfromtxt(infoFile, skip_header=1, usecols=0)
nFiles = np.genfromtxt(infoFile, skip_header=1, usecols=1)

print "      Found " + str(nFiles) + " files."
print "      Found " + str(nTetrads) + " tetrads."
print ""

# Initialize numpy arrays
data = [ structtype() for i in range(nFiles) ]

for i in range(nFiles):
  data[i].I1 = np.zeros(nTetrads)
  data[i].I2 = np.zeros(nTetrads)
  data[i].I3 = np.zeros(nTetrads)
  data[i].shape = np.zeros(nTetrads)
  data[i].var = np.zeros(nTetrads)
  data[i].RoG = np.zeros(nTetrads)

time = np.zeros(nFiles)

# Enumerate the columns
RoGCol = (0)
varCol = (1)
shapeCol = (2)
maxEig = (3)
medEig = (4)
minEig = (5)

# Loop over all output files, pull data to structures
for ff, fname in enumerate(files):
  # Pull time from filename
  ftime = fname.split("/")[-1]
  time[ff] = float(ftime[9:])

  # Pull data
  data[ff].RoG = np.genfromtxt(fname, skip_header=1, usecols=RoGCol)
  data[ff].var = np.genfromtxt(fname, skip_header=1, usecols=varCol)
  data[ff].shape = np.genfromtxt(fname, skip_header=1, usecols=shapeCol)

  data[ff].I1 = np.genfromtxt(fname, skip_header=1, usecols=maxEig)
  data[ff].I2 = np.genfromtxt(fname, skip_header=1, usecols=medEig)
  data[ff].I3 = np.genfromtxt(fname, skip_header=1, usecols=minEig)

  if ff > 99:
    break

  txt = "      File " + str(ff) + " of " + str(nFiles)
  print txt

print ""
time = time - time[0]

# Initialize histogram bins
nBins = float(35)
weights = np.ones(nTetrads)/nTetrads

# Plot specs
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=11, numpoints=3)
#plt.rc('legend', framealpha=0.7)
plt.rc('lines', markersize=4)
labelx = -0.17

# Initial Timestep
initial = plt.figure(figsize=(5,4))
gs = gridspec.GridSpec(2,3)

axEig = initial.add_subplot(gs[0,:])

n0, bins0, tmp0 = plt.hist(data[0].I1, bins=nBins, weights=weights, 
  normed=False, color="blue", alpha=0.6)
n1, bins1, tmp1 = plt.hist(data[0].I2, bins=nBins, weights=weights, 
  normed=False, color="green", alpha=0.6)
n2, bins2, tmp2 = plt.hist(data[0].I3, bins=nBins, weights=weights, 
  normed=False, color="red", alpha=0.6)

axEig.set_xlim([0,1])
axEig.set_ylim([0,0.12])
axEig.set_xlabel(r'$I_k$')
axEig.set_ylabel('Probability')
axEig.legend([r'$I_1$', r'$I_2$', r'$I_3$'])
axEig.set_xticks(np.linspace(0, 1, 11))
axEig.xaxis.set_minor_locator(AutoMinorLocator())
axEig.yaxis.set_minor_locator(AutoMinorLocator())
axEig.tick_params(which='major', length=6)
axEig.tick_params(which='minor', length=3)

x0 = bins0[np.argmax(n0)]
y0 = np.max(n0)
axEig.plot(x0, y0, 'o')
axEig.text(x0 - 0.1, y0 + 0.01, '(%.2f' % x0 + ', %.2f' % y0 + ')', fontsize=10)

x1 = bins1[np.argmax(n1)]
y1 = np.max(n1)
axEig.plot(x1, y1, 'o')
axEig.text(x1 - 0.1, y1 + 0.005, '(%.2f' % x1 + ', %.2f' % y1 + ')', fontsize=10)

x2 = bins2[np.argmax(n2)]
y2 = np.max(n2)
axEig.plot(x2, y2, 'o')
axEig.text(x2 - 0.1, y2 + 0.015, '(%.2f' % x2 + ', %.2f' % y2 + ')', fontsize=10)

# shape
axShape = plt.subplot(gs[1,0])
plt.hist(data[0].shape, weights=weights, normed=False, color='black', alpha=0.6)

axShape.set_xlabel(r'$S$')
axShape.set_ylabel('Probability')
axShape.set_xlim([-.15, .15])
axShape.set_xticks([-0.1, 0, 0.1])

# var
axVar = plt.subplot(gs[1,1])
plt.hist(data[0].var, weights=weights, normed=False, color='black', alpha=0.6)

axVar.set_xlabel(r'$\Delta$')
axVar.set_xlim([0, 0.15])
axVar.set_xticks([0, .15])

# rog
axR = plt.subplot(gs[1,2])
plt.hist(data[0].RoG, weights=weights, normed=False, 
  color='black', alpha=0.6)

axR.set_xlabel(r'$R_g [mm]$')
axR.set_xlim([2.4, 6])
axR.set_xticks([2.5,3.5,4.5,5.5])

# SAVE
plt.tight_layout()
imgname = imgdir + "shape_initial"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Time evolution of principal axes
timePrAx = plt.figure()

nPlot = 6
color = np.linspace(1, 0.2, num=nPlot)
#tstepplot = np.linspace(0,nFiles-1,nPlot)
tstepplot = np.linspace(0,100,nPlot)
tstepplot = tstepplot.astype(int)
legText = ['']*np.size(tstepplot)

i1_ax = timePrAx.add_subplot(311)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].I1, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  i1_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

i1_ax.set_xlim([0,1])
i1_ax.set_ylim([0, 0.3])
i1_ax.set_yticks([0, 0.1, 0.2, 0.3])
i1_ax.set_ylabel(r'$P(I_1)$')
i1_ax.tick_params(axis='x', labelbottom='off')

i3_ax = timePrAx.add_subplot(313)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].I3, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  i3_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

i3_ax.set_xlim([0,1])
i3_ax.set_ylim([0, 0.3])
i3_ax.set_yticks([0, 0.1, 0.2, 0.3])
i3_ax.set_ylabel(r'$P(I_3)$')
i3_ax.set_xlabel(r'$I_k$')

i2_ax = timePrAx.add_subplot(312)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].I2, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  i2_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)
  legText[i] = str(time[nn])

i2_ax.set_xlim([0,1])
i2_ax.set_ylim([0, 0.3])
i2_ax.set_yticks([0, 0.1, 0.2, 0.3])
i2_ax.set_ylabel(r'$P(I_2)$')
i2_ax.tick_params(axis='x', labelbottom='off')


#i2_ax.legend(legText, bbox_to_anchor=(0.6, 1.), loc=2, borderaxespad=0, 
#  title='Time [ms]', mode='expand')
i2_ax.legend(legText, bbox_to_anchor=(0.625, 1.15), loc=2, title='Time [ms]')

# SAVE
imgname = imgdir + "shape_principal_evolution"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## tIME EVOLUTION OF SHAPE ##
timeShape = plt.figure(figsize=(4,4))

shape_ax = timeShape.add_subplot(313)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].shape, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  shape_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

shape_ax.set_xlim([-0.25, 2])
shape_ax.set_ylim([0, 0.5])
shape_ax.set_xticks([-0.25, 0, 0.5, 1.0, 1.5, 2.0])
shape_ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
shape_ax.set_xlabel(r'$S$')
shape_ax.set_ylabel(r'$P(S)$')

var_ax = timeShape.add_subplot(312)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].var, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  var_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)
  legText[i] = str(time[nn])

var_ax.set_xlim([0,1])
var_ax.set_ylim([0, 0.4])
var_ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
var_ax.set_xlabel(r'$\Delta$')
var_ax.set_ylabel(r'$P(\Delta)$')

rg_ax = timeShape.add_subplot(311)
for i,nn in enumerate(tstepplot):
  y,edges = np.histogram(data[nn].RoG, weights=weights, normed=False)
  centers = 0.5*(edges[1:] + edges[:-1])
  rg_ax.plot(centers, y, 'ko-', alpha=color[i], linewidth=2.0)

rg_ax.set_ylim([0, 0.3])
rg_ax.set_yticks([0, 0.1, 0.2, 0.3])
rg_ax.set_xlabel(r'$R_g [mm]$')
rg_ax.set_ylabel(r'$P(R_g)$')

var_ax.legend(legText, bbox_to_anchor=(1, 0.5), loc=10,  title='Time [ms]')

# SAVE
plt.tight_layout()
imgname = imgdir + "shape_measures_evolution"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
