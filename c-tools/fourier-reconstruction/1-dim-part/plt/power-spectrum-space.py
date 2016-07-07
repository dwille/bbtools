#!/usr/bin/env python2

from setup import *

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.fftpack as scifft

os.system("clear")

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

vfAutoCorrReal = np.zeros((nz, nt))
vfAutoCorrImag = np.zeros((nz, nt))
vfPowerSpec = np.zeros((nz,nt))
#TODO: would make sense to do this starting at everyz and averaging?
for tt, tval in enumerate(time):
  (vfAutoCorrReal[:,tt], vfAutoCorrImag[:,tt], vfPowerSpec[:,tt]) \
    = AutoCorrelationSpaceFFT(vFrac[:,tt])
  vfAutoCorrReal[:,tt] /= vfAutoCorrReal[0,tt]

  #vfPowerSpec[:,tt] = np.absolute(scifft.fft(vfAutoCorrReal[:,tt]))**2
  # TODO: aren't I already doing this work in AutoCorrelationSpaceFFT?
  # TODO: actually, aren't i already doing this in the c program...
 
vfPowerSpec /= (nz*nz)
vfPowerSpecMean = np.mean(vfPowerSpec, 1)

# Make sure imag part is negligible
if (np.max(vfAutoCorrImag) > 1e-15) or (np.min(vfAutoCorrImag) < -1e-15):
  print "Imaginary part is not numerical zero:"
  print "   Max = %.2e" % np.max(vfAutoCorrImag)
  print "   Min = %.2e" % np.min(vfAutoCorrImag)
  sys.exit()

# # Plot # #
vfFig = plt.figure(figsize=(4,8))
# Plot vFrac
vfFig.add_subplot(411)
plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]],
  cmap='bwr')

# Plot autocorrelation (vertical slices at constant time)
vfFig.add_subplot(412)
print np.max(vfAutoCorrReal)
print np.min(vfAutoCorrReal)
plt.imshow(vfAutoCorrReal, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]],
  vmin=-1., vmax=1.,
  cmap='bwr')

# Plot power spectrum as a colormap
vfFig.add_subplot(413)
plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
   extent=[time[0]-0.5*dt, time[-1]+0.5*dt, 0-0.5, np.size(freq)+0.5])
#plt.plot(vfPowerSpecMean, freq, "ko:")
plt.ylim([-1, 10])
plt.ylabel(r"$k$")


 
# Plot MEAN of all power spectrums averaged over time
vfAvgPowerSpec = np.mean(vfPowerSpec, 1)
vfFig.add_subplot(414)
plt.semilogx(vfAvgPowerSpec, np.arange(0,np.size(freq)), 'ko:')
plt.xlabel('PS')
#plt.xlim([0, 0.02])
plt.ylabel(r"$k$")
plt.ylim([0, 20])
plt.xlim([1e-15, 1e2])


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
