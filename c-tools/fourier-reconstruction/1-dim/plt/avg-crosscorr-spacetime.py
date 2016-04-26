#!/usr/bin/env python2

import sys, os
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# stackoverflow 16044491 statistical scaling of autocorrelation using numpy.fft
def CrossCorrelation(x1,x2):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  x2 = np.asarray(x2)
  y2 = x2 - x2.mean()
  result = np.correlate(y2,y1,mode="full")
  result = result[len(result)/2:]
  return result

def CrossCorrelationFFT(x1,x2):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  x2 = np.asarray(x2[::-1])
  y2 = x2 - x2.mean()
  result = signal.fftconvolve(y2,y1,mode="full")
  # Re-flip result array 
  result = result[::-1]
  result = result[len(result)/2:]
  return result

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
nDensFile = datadir + "number-density"
vFracFile = datadir + "volume-fraction"
upFile = datadir + "part-u"
vpFile = datadir + "part-v"
wpFile = datadir + "part-w"

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=1)[1:] / 1000
print "      Steady state time set to: " + str(time[ts])
time = time[ts:] - time[ts]
nt = np.size(time)

evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] # / partR
nz = np.size(evalZ)
dz = evalZ - evalZ[0]

# Find output data -- each column is a different time
#numDens = np.genfromtxt(nDensFile).T[:,ts:]
vFrac = np.genfromtxt(vFracFile).T[:,ts:]
#up = np.genfromtxt(upFile).T[:,ts:]
#vp = np.genfromtxt(vpFile).T[:,ts:]
wp = np.genfromtxt(wpFile).T[:,ts:]

## Crosscorrelation of volume fraction ##
# each of (0,nz,nt) is a correlation at a starting zs
# we will average these over (nz,:,:)
vfCrossCorr = np.zeros((nz,nz,nt))
for zs in np.arange(0,nz):
  for zz, zval in enumerate(evalZ):
    # correctly loop through domain
    if zs + zz >= nz:
      zInd = zs + zz - nz
    else:
      zInd = zs + zz

    # length of result is ceil(length(time)/2)
    vfCrossCorr[zs,zz,:] = CrossCorrelationFFT(vFrac[zs,:], vFrac[zInd,:])

  vfCrossCorr[zs,:,:] /= vfCrossCorr[zs,0,0]

vfCrossCorr = np.mean(vfCrossCorr, 0)

## Save data to file
savefile = datadir + "z-averaged-vfrac-xcorr"
with open(savefile, 'wb') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerows(vfCrossCorr)

## Crosscorrelation of wp ##
# each of (0,nz,nt) is a correlation at a starting zs
# we will average these over (nz,:,:)
wpCrossCorr = np.zeros((nz,nz,nt))
for zs in np.arange(0,nz):
  for zz, zval in enumerate(evalZ):
    # correctly loop through domain
    if zs + zz >= nz:
      zInd = zs + zz - nz
    else:
      zInd = zs + zz

    # length of result is ceil(length(time)/2)
    wpCrossCorr[zs,zz,:] = CrossCorrelationFFT(wp[zs,:],wp[zInd,:])
  wpCrossCorr[zs,:,:] /= wpCrossCorr[zs,0,0]

wpCrossCorr = np.mean(wpCrossCorr, 0)

## Save data to file
savefile = datadir + "z-averaged-wp-xcorr"
with open(savefile, 'wb') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerows(wpCrossCorr)
