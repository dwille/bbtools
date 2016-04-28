#!/usr/bin/env python2

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# stackoverflow 16044491 statistical scaling of autocorrelation using numpy.fft

def AutoCorrelation(x1):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  # c_{x1,x2}[k] = sum_n x1[n+k] * conj(x2[n])
  result = np.correlate(y1,y1,mode="full")
  result = result[len(result)/2:]
  return result

def AutoCorrelationFFT(x1):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  result = signal.fftconvolve(y1[::-1],y1,mode="full")
  # Reflip array
  result = result[::-1]
  result = result[len(result)/2:]
  result /= result[0]
  return result

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                 Autocorrelate"
print ""

# SIMULATION PARAMETERS
partR = 2.1
ts = 500

# DEVEL
root = "/home/dwille/bbtools/c-tools/fourier-reconstruction/1-dim/"
simdir = "sim/"
datadir = root + simdir + "data/reconstruct-1D/"

# MARCC
#root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
#simdir = raw_input("      Simulation directory: ")
#if not simdir.endswith('/'):
#  simdir = simdir + '/'
#datadir = root + simdir + "data/reconstruct-1D/"

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
nDensFile = datadir + "number-density"
vFracFile = datadir + "volume-fraction"
upFile = datadir + "part-u"
vpFile = datadir + "part-v"
wpFile = datadir + "part-w"

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=1)[1:] / 1000
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

# Autocorrelation of volume fraction
vfFig = plt.figure()
vfAutoCorr = np.zeros((nz, nt))
vfMaxima = np.zeros((nt,3))
vfFirstMaxima = np.zeros((nt,3))
for tt,tval in enumerate(time):

  # length of result is ceil(length(time)/2)
  vfAutoCorr[:,tt] = AutoCorrelationFFT(vFrac[:,tt])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(vfAutoCorr[:,tt]))) < 0).nonzero()[0] + 1
  maxima = vfAutoCorr[maximaLoc,tt]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfMaxima[tt,0] = time[tt]
    vfMaxima[tt,1] = np.nan
    vfMaxima[tt,2] = np.nan
  else:
    vfMaxima[tt,0] = time[tt]
    vfMaxima[tt,1] = dz[maximaLoc[maxInd]]
    vfMaxima[tt,2] = vfAutoCorr[maximaLoc[maxInd],tt]

    vfFirstMaxima[tt,0] = time[tt]
    vfFirstMaxima[tt,1] = dz[maximaLoc[0]]
    vfFirstMaxima[tt,2] = vfAutoCorr[maximaLoc[0],tt]

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()
plt.plot(vfMaxima[:,0], vfMaxima[:,1], '--', color='0.75')
plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'k--')

plt.xlabel(r"$t - t_0\ [s]$")
plt.ylabel(r"$\Delta z\ [mm]$",rotation=0)
plt.title(r"$\langle \alpha(z) \alpha(z + \Delta z) \rangle $")
plt.xlim([time[0], time[-1]])
plt.ylim([dz[0], dz[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

imgname = imgdir + "autocorr-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Autocorrelation of particle velocity - space
wpFig = plt.figure()
wpAutoCorr = np.zeros((nz,nt))
wpMaxima = np.zeros((nz,3))
wpFirstMaxima = np.zeros((nt,3))
for tt,tval in enumerate(time):

  # length of result is ceil(length(time)/2)
  wpAutoCorr[:,tt] = AutoCorrelationFFT(wp[:,tt])

  # find maxima
#  maximaLoc = (np.diff(np.sign(np.diff(wpAutoCorr[:,tt]))) < 0).nonzero()[0] + 1
#  maxima = wpAutoCorr[maximaLoc,tt]
#  maxInd = np.argmax(maxima)
#  if np.size(maximaLoc) == 0:
#    wpMaxima[tt,0] = time[tt]
#    wpMaxima[tt,1] = np.nan
#    wpMaxima[tt,2] = np.nan
#  else:
#    wpMaxima[tt,0] = time[tt]
#    wpMaxima[tt,1] = dz[maximaLoc[maxInd]]
#    wpMaxima[tt,2] = wpAutoCorr[maximaLoc[maxInd],tt]
#
#    wpFirstMaxima[tt,0] = time[tt]
#    wpFirstMaxima[tt,1] = dz[maximaLoc[0]]
#    wpFirstMaxima[tt,2] = wpAutoCorr[maximaLoc[0],tt]

plt.imshow(wpAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()

#plt.plot(wpMaxima[:,0], wpMaxima[:,1], 'ko--')

plt.xlabel(r"$t - t_0\ [s]$")
plt.ylabel(r"$\Delta z\ [mm]$",rotation=0)
plt.title(r"$\langle w_p(z) w_p(z + \Delta z) \rangle $")
plt.xlim([time[0], time[-1]])
plt.ylim([dz[0], dz[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

imgname = imgdir + "autocorr-space-wp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
