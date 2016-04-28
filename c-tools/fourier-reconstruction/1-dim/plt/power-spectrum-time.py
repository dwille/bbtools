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
print "              Power Spectrum -- Time"
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

# Find output data -- each column is a different time
vFrac = np.genfromtxt(vFracFile).T[:,ts:]


dt = np.mean(np.diff(time))
samplingRate = 1./dt
freq = scifft.fftfreq(nt, dt)
# Autocorrelation of volume fraction
vfAutoCorr = np.zeros((nz, nt))
vfPowerSpec = np.zeros((nz,nt))
for zz, zval in enumerate(evalZ):
  # length of result is ceil(length(time)/2)
  vfAutoCorr[zz,:] = AutoCorrelationFFT(vFrac[zz,:])
  vfPowerSpec[zz,:] = np.absolute(scifft.fft(vfAutoCorr[zz,:]))**2

vfPowerSpec /= (nt*nt)
vfAutoCorrMean = np.mean(vfAutoCorr, 0)
vfPowerSpecMean = np.mean(vfPowerSpec, 0)

plt.rc('figure', figsize=(4,5))
vfFig = plt.figure()
vfFig.add_subplot(311)
plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])

vfFig.add_subplot(312)
plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
  extent=[0, np.max(freq*2), evalZ[0], evalZ[-1]])
#plt.plot(freq, np.sqrt(vfPowerSpecMean), "ko:")
plt.xlim([0, 5])


#vfAutoCorr = AutoCorrelationFFT(vFrac[zs,:])

# Plot autocorrelation
#vfFig.add_subplot(211)
#plt.plot(time, vfAutoCorr)
#plt.xlabel('Tau')
#plt.ylabel('Autocorrelation')
#plt.yticks([-1, -.5, 0, 0.5, 1])
#
vfFig.add_subplot(313)
# Method 1
#fftSignal = np.abs(scifft.fft(vfAutoCorr[0,:]))**2
#W = scifft.fftfreq(nt, dt)
#plt.plot(W, np.sqrt(fftSignal), "ko:")

zval = 480
# Method 2
freqs = np.fft.fftfreq(nt, dt)
idx = np.argsort(freqs)

ps = np.absolute(np.fft.fft(vfAutoCorr[0,:]))**2/(nt*nt)
plt.plot(freqs[idx], ps[idx], "ro:")

ps = np.absolute(np.fft.fft(vfAutoCorr[zval,:]))**2/(nt*nt)
plt.plot(freqs[idx], ps[idx], "bo:")

plt.legend(["z = %.2f" % evalZ[0], "z = %2.f" % evalZ[zval]])


plt.xlabel("Freq")
plt.ylabel(r"$PS$")

plt.xlim([0, 5])
plt.ylim([0, 0.005])

#vfFig.add_subplot(313)
#f, Pxx_den = signal.periodogram(vfAutoCorr, 1./dt)
#plt.bar(f, Pxx_den)
#plt.xlim([0,5])
#plt.ylabel("log(sqrt(FT(autocorrelation)))")

imgname = imgdir + "power-spectrum-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
