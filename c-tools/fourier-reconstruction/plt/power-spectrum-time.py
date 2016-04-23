#!/usr/bin/env python2

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.fftpack import rfft, irfft, fftfreq

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
time = np.genfromtxt(infoFile, skip_footer=1)[1:]
time = time[ts:] - time[ts]
nt = np.size(time)

evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] / partR
nz = np.size(evalZ)

# Find output data -- each column is a different time
#numDens = np.genfromtxt(nDensFile).T[:,ts:]
vFrac = np.genfromtxt(vFracFile).T[:,ts:]
#up = np.genfromtxt(upFile).T[:,ts:]
#vp = np.genfromtxt(vpFile).T[:,ts:]
#wp = np.genfromtxt(wpFile).T[:,ts:]

# Plot specs
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
#plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,5))
plt.rc('legend', fontsize=11, numpoints=3)
plt.rc('lines', markersize=2)
plt.rc('savefig', dpi=250)
labelx = -0.17

# Autocorrelation of volume fraction
vfFig = plt.figure()
#vfAutoCorr = np.zeros((nz, nt))
#vfPowerSpec = np.zeros((nz,3))
#for zz, zval in enumerate(evalZ):
#
#  # length of result is ceil(length(time)/2)
#  vfAutoCorr[zz,:] = AutoCorrelationFFT(vFrac[zz,:])
#  vfPowerSpec[zz,:] = np.abs(np.fft.fft(vfAutoCorr[zz,:]))**2
vfAutoCorr = AutoCorrelationFFT(vFrac[0,:])
dt = np.mean(np.diff(time))/1000
W = fftfreq(nt, dt)
fftSignal = rfft(vfAutoCorr)

vfFig.add_subplot(211)
plt.plot(time, vfAutoCorr)

#vfFig.add_subplot(312)
#plt.bar(W, fftSignal, 0.25)
#plt.xlim([0,5])

vfFig.add_subplot(212)
plt.bar(W, np.log(np.sqrt(fftSignal**2)), 0.25)
plt.yscale('log')
plt.xlim([0,5])
plt.ylabel("log(sqrt(FT(autocorrelation)))")

#vfFig.add_subplot(313)
#f, Pxx_den = signal.periodogram(vfAutoCorr, 1./dt)
#plt.bar(f, Pxx_den)
#plt.xlim([0,5])
#plt.ylabel("log(sqrt(FT(autocorrelation)))")

imgname = imgdir + "power-spectrum-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
