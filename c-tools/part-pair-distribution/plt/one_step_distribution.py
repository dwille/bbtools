#!/usr/bin/env python2

import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt

from floating_polar import fractional_polar_axes

print ""
print " ---- Particle Pair Distribution Plotting Utility ---- "
print ""

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

# SIMULATION PARAMETERS
partR = 2.1
ts = 500

# DEVEL
root = "/home/dwille/bbtools/c-tools/part-pair-distribution/"
simdir = "sim/"
datadir = root + simdir + "data-part-pair/"

# MARCC
#root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
#simdir = raw_input("      Simulation directory: ")
#if not simdir.endswith('/'):
#  simdir = simdir + '/'
#datadir = root + simdir + "data-reconstruct/"

print "      Sim root directory set to: " + root

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(datadir):
  print "      " + datadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + simdir + "/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# set up output file paths
infoFile = datadir + "info"

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=2)[1:]
#time = time[ts:] - time[ts]
nt = np.size(time)

evalR = np.genfromtxt(infoFile, skip_header=1,skip_footer=1)[1:] / partR
nr = np.size(evalR)

evalTh = np.genfromtxt(infoFile, skip_header=2)[1:]
nth = np.size(evalTh)

# Find output data -- each column is a different time
files = sorted_nicely(glob.glob(datadir + "part-pair-*"))
nFiles = len(files)

print "      Found " + str(nFiles) + " files."

# Loop and pull
data = np.zeros((nt, nr, nth));
for ff, fname in enumerate(files):
  data[ff,:,:] = np.genfromtxt(fname)

# Plot specs
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
#plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=11, numpoints=3)
plt.rc('lines', markersize=2)
plt.rc('savefig', dpi=250)
labelx = -0.17

## Plot ###
theta,rad = np.meshgrid(evalTh*180/np.pi - 90, evalR)
X = theta
Y = rad

f1 = plt.figure()
a1 = fractional_polar_axes(f1, thlim=(0,90), step=(15,0.2),
  theta_offset=0)
  
a1.pcolormesh(X,Y,data[0,:,:])
#ax = fig.add_subplot(111, projection="polar")

# TODO: clip polar plot
#p1 = ax.pcolormesh(X,Y,data[0,:,:])
#plt.colorbar(p1)
#ax.set_theta_direction(-1)
#ax.set_theta_zero_location("N")


imgname = imgdir + "part-pair-disti"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
