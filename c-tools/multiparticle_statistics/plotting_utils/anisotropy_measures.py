#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import re

def sorted_nicely( l ):
    """ Sorts the given iterable in a natural way

    Required arguments:
    l -- The iterable to be sroted.

    courtesy stackoverflow/2669059
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

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

# Initialize lists
varList = list()
detList = list()
R2List = list()

# Sort output files and find number within time range
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

print "Found " + str(n) + " files in time range"
print ""


nFiles = n
# Initialize mean/std arrays
varMean = np.zeros(nFiles)
varStd = np.zeros(nFiles)
shapeMean = np.zeros(nFiles)
shapeStd = np.zeros(nFiles)
RMean = np.zeros(nFiles)
RStd = np.zeros(nFiles)
time = np.zeros(nFiles)

# Initialize histogram bins
nBins = float(100)
rMax = 20
shapeHist = np.zeros([nBins, nFiles])
varHist = np.zeros([nBins, nFiles])
rHist = np.zeros([nBins, nFiles])

shapeEdges = np.linspace(-0.25, 2.0, nBins + 1)
varEdges = np.linspace(0, 1, nBins + 1)
rEdges = np.linspace(0, rMax, nBins + 1)

shapeCenters = shapeEdges[0:nBins] + 0.5*2.25/nBins
varCenters = varEdges[0:nBins] + 0.5*1/nBins
rCenters = rEdges[0:nBins] + 0.5*rMax/nBins

# Loop over all output files
i = 0;
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('tetrad-data-'):
    ftime = ftime[12:]

  if float(ftime) >= float(ts) and float(ftime) <= float(te):
    time[i] = float(ftime)
    ifile = open(fname)
    ifile.readline()

    # Loop over all lines in files
    for line in ifile:
      buf = line.split()
      detList.append(float(buf[0]))
      varList.append(float(buf[1]))
      R2List.append(float(buf[2]))

    ifile.close()

    # Pull values to np arrays
    det = np.array(detList)
    var = np.array(varList)
    R2 = np.array(R2List)

    # Calculate shape and variance measure for each tetrad
    tetradShape = 27*det/(R2*R2*R2)
    tetradEVar = 1.5*var/(R2*R2)

    # Bin the data
    shapeHist[:,i], tmp = np.histogram(tetradShape, bins=shapeEdges,
      range=[-0.25, 2], density=True)
    varHist[:,i], tmp = np.histogram(tetradEVar, bins=varEdges,
      range=[-0.25, 2], density=True)
    rHist[:,i], tmp = np.histogram(np.sqrt(R2), bins=rEdges, density=True)

    # At each tstep, average measures over all tetrads
    shapeMean[i] = np.mean(tetradShape)
    shapeStd[i] = np.std(tetradShape)

    varMean[i] = np.mean(tetradEVar)
    varStd[i] = np.std(tetradEVar)

    RMean[i] = np.mean(np.sqrt(R2))
    RStd[i] = np.std(np.sqrt(R2))

    i += 1


# Plot over time
measures = plt.figure()
measures.suptitle('Anisotropy Measures', fontsize=20)
labelx = -0.10

# Radius of Gyration
rg_ax = measures.add_subplot(311)
rg_ax.plot(time, RMean, 'wo-', linewidth=1.5)
rg_ax.plot(time, RMean + RStd, 'k--', color='0', linewidth=2.5)
rg_ax.plot(time, RMean - RStd, 'k--', color='0', linewidth=2.5)

rg_ax.set_ylabel("Radius of Gyration")
rg_ax.set_ylim([0,16])
rg_ax.yaxis.set_label_coords(labelx, 0.5)
rg_ax.set_xticks(np.linspace(time[0],time[nFiles-1],6))
rg_ax.set_yticks(np.linspace(0,16,5))

# Shape Factor
sf_ax = measures.add_subplot(312)
sf_ax.plot(time, shapeMean, 'wo-')
sf_ax.plot(time, shapeMean + shapeStd, 'k--', color='0', linewidth=2.5)
sf_ax.plot(time, shapeMean - shapeStd, 'k--', color='0', linewidth=2.5)

sf_ax.set_ylabel("Shape Factor")
sf_ax.yaxis.set_label_coords(labelx, 0.5)
sf_ax.set_xlim([time[0], time[nFiles-1]])
sf_ax.set_ylim([-0.25,2])
sf_ax.set_xticks(np.linspace(time[0],time[nFiles-1],6))
sf_ax.set_yticks([-0.25,0,0.5,1.0,1.5,2.0])

# Variance
var_ax = measures.add_subplot(313)
var_ax.plot(time, varMean, 'wo-')
var_ax.plot(time, varMean + varStd, 'k--', color='0', linewidth=2.5)
var_ax.plot(time, varMean - varStd, 'k--', color='0', linewidth=2.5)
 
var_ax.set_xlabel("Time [ms]")
var_ax.set_ylabel("Eigenvalue Variance")
var_ax.yaxis.set_label_coords(labelx, 0.5)
var_ax.set_xlim([time[0], time[nFiles-1]])
var_ax.set_ylim([0,1])
sf_ax.set_xticks(np.linspace(time[0],time[nFiles-1],6))
var_ax.set_yticks(np.arange(0,1.1,0.2))

# Surface plots

# Radius of Gyration
rImg = rg_ax.imshow(rHist, origin='lower', aspect='auto', interpolation='none',
         extent=[np.amin(time), np.amax(time), 0, rMax], cmap=cm.jet)
pos1 = rg_ax.get_position()
pos2 = [pos1.x0, pos1.y0 - 0.01, pos1.width*0.9, pos1.height]
rg_ax.set_position(pos2)

# Shape Factor
sImg = sf_ax.imshow(shapeHist, origin='lower', aspect='auto', 
  interpolation='none', extent=[np.amin(time), np.amax(time), -0.25, 2],
  cmap=cm.jet)
pos1 = sf_ax.get_position()
pos2 = [pos1.x0, pos1.y0 - 0.01, pos1.width*0.9, pos1.height]
sf_ax.set_position(pos2)

# Eig Variance
var_ax = measures.add_subplot(313)
vImg = var_ax.imshow(varHist, origin='lower', aspect='auto', 
  interpolation='none', extent=[np.amin(time), np.amax(time), 0, 1], 
  cmap=cm.jet)
pos1 = var_ax.get_position()
pos2 = [pos1.x0, pos1.y0 - 0.01, pos1.width*0.9, pos1.height]
var_ax.set_position(pos2)

# Colorbar
cbaxes = measures.add_axes([0.9, 0.1 - 0.02, 0.03, 0.8])
cbar = measures.colorbar(rImg, cax = cbaxes)
cbar.ax.set_title('Probability')

#plt.tight_layout()
plt.show()

# Shape factor
# # Set up Latex plot rendering
# fs = 12
# #plt.rc('text', usetex=True)
# #plt.rc('font', family='serif')
# #plt.rc('legend', fontsize=fs)
# #xticklabels = plt.getp(plt.gca(), 'xticklabels')
# #yticklabels = plt.getp(plt.gca(), 'yticklabels')
# #plt.setp(xticklabels, fontsize=fs)
# #plt.setp(yticklabels, fontsize=fs)
