#!/usr/bin/env python2

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.fftpack as scifft

os.system("clear")

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
print "             Power Spectrum -- Space"
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
vFracFile = datadir + "volume-fraction"

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=1)[1:] / 1000
time = time[ts:] - time[ts]
nt = np.size(time)

evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] # / partR
nz = np.size(evalZ)
dz = evalZ - evalZ[0]

# Find output data -- each column is a different time
vFrac = np.genfromtxt(vFracFile).T[:,ts:]

freq = scifft.fftfreq(nz, np.mean(np.diff(dz)))
# Autocorrelation of volume fraction
vfAutoCorr = np.zeros((nz, nt))
vfPowerSpec = np.zeros((nz,nt))
for tt, tval in enumerate(time):
  # length of result is ceil(length(time)/2)
  vfAutoCorr[:,tt] = AutoCorrelationFFT(vFrac[:,tt])
  vfPowerSpec[:,tt] = np.abs(scifft.fft(vfAutoCorr[:,tt]))**2

vfPowerSpec /= (nz*nz)
vfAutoCorrMean = np.mean(vfAutoCorr, 1)
vfPowerSpecMean = np.mean(vfPowerSpec, 1)

plt.rc('figure', figsize=(4,5))
vfFig = plt.figure()
vfFig.add_subplot(311)
plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])

vfFig.add_subplot(312)
plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], 0, np.max(freq)])
#plt.plot(vfPowerSpecMean, freq, "ko:")
plt.ylim([0, 0.05])


#vfAutoCorr = AutoCorrelationFFT(vFrac[zs,:])

# Plot autocorrelation
#vfFig.add_subplot(211)
#plt.plot(time, vfAutoCorr)
#plt.xlabel('Tau')
#plt.ylabel('Autocorrelation')
#plt.yticks([-1, -.5, 0, 0.5, 1])
#
vfFig.add_subplot(313)
## Method 1
#fftSignal = np.abs(scifft.fft(vfAutoCorr[0,:]))**2
#W = scifft.fftfreq(nt, dt)
#plt.plot(W, np.sqrt(fftSignal), "ko:")
#
tval = 1200
# Method 2
freqs = np.fft.fftfreq(nz, np.mean(np.diff(dz)))
idx = np.argsort(freqs)

ps = np.abs(np.fft.fft(vfAutoCorr[:,0]))**2/(nz*nz)
plt.plot(ps[idx], freqs[idx], "ro:")

ps = np.abs(np.fft.fft(vfAutoCorr[:,tval]))**2/(nz*nz)
plt.plot(ps[idx], freqs[idx], "bo:")

plt.legend(["t = %.2f" % time[0], "t = %.2f" % time[tval]],
  loc="lower right")

plt.xlabel("PS")
plt.ylabel("lambda?")
#
plt.xlim([0, 0.02])
plt.ylim([0, 0.05])

#vfFig.add_subplot(313)
#f, Pxx_den = signal.periodogram(vfAutoCorr, 1./dt)
#plt.bar(f, Pxx_den)
#plt.xlim([0,5])
#plt.ylabel("log(sqrt(FT(autocorrelation)))")

imgname = imgdir + "power-spectrum-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
