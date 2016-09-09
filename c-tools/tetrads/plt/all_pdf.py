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

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line

print "Anisotropy Measures Plotting Utility"
print "   Plot measures for all simulations"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print " Sim root directory set to: " + root

# Parameter sweep
partList = ['500', '1000', '1500', '2000']
densList = ['rho2.0', 'rho3.3', 'rho4.0', 'rho5.0']
legendText = ['']*16

# Initialize mean/std arrays
data = [ structtype() for i in range(16) ]
 
# Loop over all directory's
for pp, part in enumerate(partList):
  for dd, dens in enumerate(densList):
    # Set up path to data
    stride = 4*pp + dd
    caseDir = part + '/' + dens
    statFile = root + caseDir + '/' + 'data-tetrads/stat.dat'
    nodeFile = root + caseDir + '/' + 'data-tetrads/regularNodes'
    rawFiles = root + caseDir + '/' + 'data-tetrads/raw-data-*'

    # Find number of tetrads
    data[stride].nTetrads = file_len(nodeFile) 

    # Find raw data files and number of timesteps
    rawFiles = sorted_nicely(glob.glob(rawFiles))
    data[stride].ntSteps = len(rawFiles)

    # Pull time
    time = np.genfromtxt(statFile, skip_header=1, usecols=0)
    data[stride].time = np.zeros(np.size(time))
    data[stride].time = time - time[0]

    # Init variables
    tempG1 = np.zeros((data[stride].nTetrads, 3))
    tempG2 = np.zeros((data[stride].nTetrads, 3))
    tempG3 = np.zeros((data[stride].nTetrads, 3))

    weights = np.ones(data[stride].nTetrads)/data[stride].nTetrads
    nbins = 25

    data[stride].count_g1z = np.zeros((nbins, data[stride].ntSteps))
    data[stride].count_g2z = np.zeros((nbins, data[stride].ntSteps))
    data[stride].count_g3z = np.zeros((nbins, data[stride].ntSteps))
    data[stride].center_g1z = np.zeros((nbins, data[stride].ntSteps))
    data[stride].center_g2z = np.zeros((nbins, data[stride].ntSteps))
    data[stride].center_g3z = np.zeros((nbins, data[stride].ntSteps))

    # Loop over all files in the directory
    for ff, fname in enumerate(rawFiles):
      # Pull data to tmp variables
      tempGV1 = np.genfromtxt(fname, skip_header=1, usecols=(6,7,8))
      tempGV2 = np.genfromtxt(fname, skip_header=1, usecols=(9,10,11))
      tempGV3 = np.genfromtxt(fname, skip_header=1, usecols=(12,13,14))

      # Dot with z = (0,0,1)
      g1z = tempGV1[:,2]
      g2z = tempGV2[:,2]
      g3z = tempGV3[:,2]

      y_g1z, edges_g1z = np.histogram(g1z, bins=nbins, weights=weights, 
        normed=False)
      y_g2z, edges_g2z = np.histogram(g2z, bins=nbins, weights=weights, 
        normed=False)
      y_g3z, edges_g3z = np.histogram(g3z, bins=nbins, weights=weights, 
        normed=False)

      data[stride].count_g1z[:,ff] = y_g1z
      data[stride].count_g2z[:,ff] = y_g2z
      data[stride].count_g3z[:,ff] = y_g3z
      data[stride].center_g3z[:,ff] = 0.5*(edges_g1z[1:] + edges_g1z[:-1])
      data[stride].center_g2z[:,ff] = 0.5*(edges_g2z[1:] + edges_g2z[:-1])
      data[stride].center_g3z[:,ff] = 0.5*(edges_g3z[1:] + edges_g3z[:-1])




#    data[stride].g1z = np.zeros(np.size(time))
#    data[stride].g2z = np.zeros(np.size(time))
#    data[stride].g3z = np.zeros(np.size(time))
#
#    data[stride].g1z = np.genfromtxt(alignFile, skip_header=1, usecols=10)
#    data[stride].g2z = np.genfromtxt(alignFile, skip_header=1, usecols=11)
#    data[stride].g3z = np.genfromtxt(alignFile, skip_header=1, usecols=12)
#
#    legendText[stride] = caseDir + ': ' + str(data[stride].nnodes)

plt.rc('font', family='serif')
color = ['r', 'g', 'b', 'k']
shade = [0.4, 0.57, 0.74, 0.9]


### Shape and Strain Alignment ##
#
############
### gi_s1 ##
############
#gi_s1_Fig = plt.figure(figsize=(12,8))
#gi_s1_Fig.suptitle(r'$\langle (g_i, s_1) \rangle$', fontsize=16)
#
## n = 500
#gi_s1_500_ax = gi_s1_Fig.add_subplot(411)
#for dd in range(4):
#  i = dd
#  gi_s1_500_ax.plot(data[i].time,data[i].g1s1, linewidth=3,
#    color='k', alpha=shade[dd])
#  gi_s1_500_ax.plot(data[i].time,data[i].g2s1, linewidth=3,
#    color='r', alpha=shade[dd])
#  gi_s1_500_ax.plot(data[i].time,data[i].g3s1, linewidth=3,
#    color='b', alpha=shade[dd])
#    
#gi_s1_500_ax.set_ylabel('n = 500')
#gi_s1_500_ax.set_xlim([0, 300])
#gi_s1_500_ax.set_ylim([0, 0.7])
#gi_s1_500_ax.grid(True)
#
#plt.show()
