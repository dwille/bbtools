#!/usr/bin/env python2

import sys, os
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                 Crosscorrelate"
print ""

# SIMULATION PARAMETERS
partR = 2.1
ts = 500

# DEVEL
root = "/home/dwille/bbtools/c-tools/fourier-reconstruction/"
simdir = "sim/"
datadir = root + simdir + "data-reconstruct/"

# MARCC
#root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
#simdir = raw_input("      Simulation directory: ")
#if not simdir.endswith('/'):
#  simdir = simdir + '/'
#datadir = root + simdir + "data-reconstruct/"

print "      Sim root directory set to: " + root
print "      Particle Radius set to: " + str(partR)
print "      Steady state index set to: " + str(ts)

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(datadir):
  print "      " + datadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + simdir + "/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# set up data file paths
infoFile = datadir + "info"
vfSaveFile = datadir + "z-averaged-vfrac-xcorr"

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=1)[1:] / 1000
time = time[ts:] - time[ts]
nt = np.size(time)

evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] #/ partR
nz = np.size(evalZ)
dz = evalZ - evalZ[0]

# Pull output data -- each column is a different time
vfCrossCorr = np.genfromtxt(vfSaveFile)

# Find maxima
vfFirstMaxima = np.zeros((nz,3))
for zz in np.arange(0,nz):
  maximaLoc = (np.diff(np.sign(np.diff(vfCrossCorr[zz,:]))) < 0).nonzero()[0] + 1
  maxima = vfCrossCorr[zz,maximaLoc]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfFirstMaxima[zz,0] = np.nan
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = np.nan
  else:
    vfFirstMaxima[zz,0] = time[maximaLoc[0]]
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[0]]

# Find slope of maxima -- wavespeed
for zz in np.arange(1,nz):
  # if the time changes drastically, cut the array
  if np.abs(vfFirstMaxima[zz,0] - vfFirstMaxima[zz-1,0]) > .100:
    firstWave = vfFirstMaxima[0:zz-1,:]

tau = firstWave[:,0]
tau = tau[:,np.newaxis]
dzMax = firstWave[:,1]
dzMax = dzMax[:,np.newaxis]

# Make sure that x,y[0] = 0
tau[0] = 0
dzMax[0] = 0
x = firstWave[:,0]

# Fit curve, assume 0 intercept
p, _, _, _ = np.linalg.lstsq(tau, dzMax)
xFit = firstWave[:,0]
yFit = p[0,0]*xFit

## PLOTTING ##
# Volume Fraction
vfFig = plt.figure()
plt.imshow(vfCrossCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()
plt.plot(firstWave[:,0], firstWave[:,1], 'ko', markevery=25, 
  markerfacecolor="None")
plt.plot(xFit, yFit, '--')
cTxtString = r"$dz = %.4f\tau$" % p[0,0]
plt.text(xFit[-1], yFit[-1], cTxtString, fontsize=12)

plt.xlabel(r"$\tau\ [s]$")
plt.xlim([0, time[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

plt.ylabel(r"$dz\ [mm]$", rotation=0)
plt.ylim([dz[0], dz[-1]])


plt.title(r"$\langle \alpha(t,z) \alpha(t + \tau, z + \Delta z) \rangle$")

imgname = imgdir + "avg-crosscorr-spacetime-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# # Particle vertical velocity
# wpFig = plt.figure()
# plt.imshow(wpCrossCorr, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], dz[0], dz[-1]])
# plt.colorbar()
# 
# plt.xlabel(r"$\tau\ [s]$")
# plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
# plt.ylabel(r"$dz\ [mm]$", rotation=0)
# plt.title(r"$\langle w_p(t,z) w_p(t + \tau, z + \Delta z) \rangle$")
# 
# imgname = imgdir + "avg-crosscorr-spacetime-wp"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
