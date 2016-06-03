#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import sys, os
os.system('clear')

## GET INFO
print " ---- Phase-Averaged Fluid Velocity Plotting Utility ---- "
print ""

# DEVEL
#root = "/home/dwille/bbtools/c-tools/phase-averaged-velocity/"
#simdir = "sim/"
#datadir = root + simdir + "data-reconstruct/"

# MARCC
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
simdir = raw_input("      Simulation directory: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'
datadir = root + simdir + "phasevel/data/"

print "      Sim root directory set to: " + root
print "      Using location: " + datadir

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(datadir):
  print "      " + datadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + simdir + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Set up output file paths
dataFile = datadir + "phaseAveragedVel"

time = np.genfromtxt(dataFile, skip_header=1, usecols=0)
uf = np.genfromtxt(dataFile, skip_header=1, usecols=1)
vf = np.genfromtxt(dataFile, skip_header=1, usecols=2)
wf= np.genfromtxt(dataFile, skip_header=1, usecols=3)

# Plot
phaseVel = plt.figure(figsize=(12,8))
phaseVel.suptitle('Phase Averaged Fluid Velocity', fontsize=20)

maxY = 1.1*np.max(wf)
minY = np.min(wf) - 0.1*maxY
ind = np.arange(0,np.size(time), 1)

# u
uAx = phaseVel.add_subplot(311)
uAx.plot(ind, uf, 'k-')
uAx.set_ylim([minY, maxY])
uAx.set_ylabel('u [m/s]')

# v
vAx = phaseVel.add_subplot(312)
vAx.plot(ind, vf, 'k-')
vAx.set_ylim([minY, maxY])
vAx.set_ylabel('v [m/s]')

# w
wAx = phaseVel.add_subplot(313)
wAx.plot(ind, wf, 'k-')
wAx.set_ylim([minY, maxY])
wAx.set_ylabel('w [m/s]')
wAx.set_xlabel('Index')

majorTicks = np.arange(0,np.size(time), 100)
minorTicks = np.arange(0,np.size(time), 50)
wAx.set_xticks(majorTicks, minor=False)
wAx.set_xticks(minorTicks, minor=True)
wAx.xaxis.grid(True, which='major')
wAx.xaxis.grid(True, which='minor')

plt.show()

indW = raw_input('Choose an index where you want to take a mean from: ')
indW = int(indW)
meanW = np.mean(wf[indW:])

print "Steady state value of wf = %.5f" % meanW
