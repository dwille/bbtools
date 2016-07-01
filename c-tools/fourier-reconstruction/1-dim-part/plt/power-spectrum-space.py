#!/usr/bin/env python2

from setup import *

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

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Autocorrelation of volume fraction
dt = np.mean(np.diff(time))
dz = evalZ - evalZ[0]
deltaZ = np.mean(np.diff(dz))
freq = scifft.fftfreq(nz, deltaZ)

vfAutoCorr = np.zeros((nz, nt))
vfPowerSpec = np.zeros((nz,nt))
# for tt, tval in enumerate(time):
#   # length of result is ceil(length(time)/2)
#   vfAutoCorr[:,tt] = AutoCorrelationFFT(vFrac[:,tt])
#   vfPowerSpec[:,tt] = np.absolute(scifft.fft(vfAutoCorr[:,tt]))**2
# 
# vfPowerSpec /= (nz*nz)
# vfAutoCorrMean = np.mean(vfAutoCorr, 1)
# vfPowerSpecMean = np.mean(vfPowerSpec, 1)
# 
# # Plot
vfFig = plt.figure()
# 
# # Plot autocorrelation (slices at constant time)
# vfFig.add_subplot(311)
# plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], dz[0], dz[-1]])

vfFig.add_subplot(211)
array = vFrac[:,0] - np.mean(vFrac[:,0])
length = np.size(array)
for i in np.arange(0,length):
  array[i] = np.sin(2.*np.pi*i/length)
plt.plot(array, np.arange(0,length), 'b')
#plt.ylim([0, 120])
#plt.plot(test, dz, 'k')


vfFig.add_subplot(212)
R = np.fft.ifft(np.fft.fft(array)*np.fft.fft(array[::-1]))
plt.plot(R, np.arange(0, np.size(R)), 'k--')

# # Plot ALL power spectrums as color map
# vfFig.add_subplot(312)
# plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
#   #extent=[time[0], time[-1], 0, np.max(freq)])
#   extent=[time[0]-0.5*dt, time[-1]+0.5*dt, 0-0.5, np.size(freq)+0.5])
# #plt.plot(vfPowerSpecMean, freq, "ko:")
# plt.ylim([-1, 10])
# plt.ylabel(r"$k$")
# 
# # Plot MEAN of all power spectrums averaged over time
# vfAvgPowerSpec = np.mean(vfPowerSpec, 1)
# 
# vfFig.add_subplot(313)
# plt.semilogx(vfAvgPowerSpec, np.arange(0,np.size(freq)), 'ko:')
# plt.xlabel('PS')
# #plt.xlim([0, 0.02])
# plt.ylabel(r"$k$")
# plt.ylim([0, 20])
# plt.xlim([1e-4, 1e-1])


#  #vfAutoCorr = AutoCorrelationFFT(vFrac[zs,:])
#  
#  # Plot autocorrelation
#  #vfFig.add_subplot(211)
#  #plt.plot(time, vfAutoCorr)
#  #plt.xlabel('Tau')
#  #plt.ylabel('Autocorrelation')
#  #plt.yticks([-1, -.5, 0, 0.5, 1])
#  #
#  vfFig.add_subplot(313)
#  ## Method 1
#  #fftSignal = np.abs(scifft.fft(vfAutoCorr[0,:]))**2
#  #W = scifft.fftfreq(nt, dt)
#  #plt.plot(W, np.sqrt(fftSignal), "ko:")
#  #
#  tval = 1200
#  # Method 2
#  freqs = np.fft.fftfreq(nz, np.mean(np.diff(dz)))
#  idx = np.argsort(freqs)
#  
#  ps = np.absolute(np.fft.fft(vfAutoCorr[:,0]))**2/(nz*nz)
#  plt.plot(ps[idx], freqs[idx], "ro:")
#  
#  ps = np.absolute(np.fft.fft(vfAutoCorr[:,tval]))**2/(nz*nz)
#  plt.plot(ps[idx], freqs[idx], "bo:")
#  
#  plt.legend(["t = %.2f" % time[0], "t = %.2f" % time[tval]],
#    loc="lower right")
#  
#  plt.xlabel("PS")
#  plt.ylabel("lambda?")
#  #
#  plt.xlim([0, 0.02])
#  plt.ylim([0, 0.05])
#  
#  #vfFig.add_subplot(313)
#  #f, Pxx_den = signal.periodogram(vfAutoCorr, 1./dt)
#  #plt.bar(f, Pxx_den)
#  #plt.xlim([0,5])
#  #plt.ylabel("log(sqrt(FT(autocorrelation)))")

imgname = imgdir + "power-spectrum-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
