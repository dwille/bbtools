#!/usr/bin/env python2

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

os.system('clear')

# stackoverflow 16044491 statistical scaling of autocorrelation using numpy.fft
def AutoCorrelation(x):
  x = np.asarray(x)
  y = x - x.mean()
  result = np.correlate(y,y,mode="full")
  result = result[len(result)/2:]
  result /= result[0]
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

evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] #/ partR
nz = np.size(evalZ)

# Find output data -- each column is a different time
numDens = np.genfromtxt(nDensFile).T[:,ts:]
vFrac = np.genfromtxt(vFracFile).T[:,ts:]
up = np.genfromtxt(upFile).T[:,ts:]
vp = np.genfromtxt(vpFile).T[:,ts:]
wp = np.genfromtxt(wpFile).T[:,ts:]

# Plot specs
#plt.rc('xtick', labelsize=10)
#plt.rc('ytick', labelsize=10)
#plt.rc('axes', labelsize=11)
##plt.rc('figure', titlesize=14)
#plt.rc('figure', figsize=(4,3))
#plt.rc('legend', fontsize=11, numpoints=3)
#plt.rc('lines', markersize=2)
#plt.rc('savefig', dpi=250)
labelx = -0.17

# Autocorrelation of volume fraction
vfFig = plt.figure()
vfAutoCorr = np.zeros((nz, nt))
vfMaxima = np.zeros((nz,3))
vfFirstMaxima = np.zeros((nz,3))
for zz, zval in enumerate(evalZ):

  # length of result is ceil(length(time)/2)
  vfAutoCorr[zz,:] = AutoCorrelationFFT(vFrac[zz,:])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(vfAutoCorr[zz,:]))) < 0).nonzero()[0] + 1
  maxima = vfAutoCorr[zz,maximaLoc]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfMaxima[zz,0] = np.nan
    vfMaxima[zz,1] = evalZ[zz]
    vfMaxima[zz,2] = np.nan
  else:
    vfMaxima[zz,0] = time[maximaLoc[maxInd]]
    vfMaxima[zz,1] = evalZ[zz]
    vfMaxima[zz,2] = vfAutoCorr[zz,maximaLoc[maxInd]]

    vfFirstMaxima[zz,0] = time[maximaLoc[0]]
    vfFirstMaxima[zz,1] = evalZ[zz]
    vfFirstMaxima[zz,2] = vfAutoCorr[zz,maximaLoc[0]]

tauVal = np.mean(vfFirstMaxima[:,0])
tauMean = np.mean(vfFirstMaxima[:,2])
print ""
print "      Mean = %.2f at tau = %.4f" % (tauMean, tauVal)
print "      omega = 2pi/tau = %.4f" % (2.*np.pi/tauVal)
print "      Freq = 1/tau = %.4f" % (1./tauVal)

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.colorbar()
#plt.plot(vfMaxima[:,0], vfMaxima[:,1], '--', color='0.75')
plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'k--')

plt.plot([tauVal, tauVal], [evalZ[0], evalZ[-1]], 'k:')
txtString = r"$\bar \tau = %.4f$" % tauVal
plt.text(1, 0, txtString, fontsize=12)

plt.xlabel(r"$\tau\ [s]$")
plt.ylabel(r"$z\ [mm]$",rotation=0)
plt.title(r"$\langle \alpha(t) \alpha(t + \tau) \rangle $")
plt.xlim([time[0], time[-1]])
plt.ylim([evalZ[0], evalZ[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

imgname = imgdir + "autocorr-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Autocorrelation of wp -- time
wpFig = plt.figure()
wpAutoCorr = np.zeros((nz, nt))
wpMaxima = np.zeros((nz,3))
for zz, zval in enumerate(evalZ):
  # length of result is ceil(length(time)/2)
  wpAutoCorr[zz,:] = AutoCorrelationFFT(wp[zz,:])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(wpAutoCorr[zz,:]))) < 0).nonzero()[0] + 1
  if np.size(maximaLoc) == 0:
    wpMaxima[zz,0] = np.nan
    wpMaxima[zz,1] = evalZ[zz]
    wpMaxima[zz,2] = np.nan
  else:
    wpMaxima[zz,0] = time[maximaLoc[0]]
    wpMaxima[zz,1] = evalZ[zz]
    wpMaxima[zz,2] = wpAutoCorr[zz,maximaLoc[0]]

plt.imshow(wpAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.colorbar()
plt.plot(wpMaxima[:,0], wpMaxima[:,1], 'k--')

plt.xlabel(r"$\tau\ [s]$")
plt.ylabel(r"$z\ [mm]$",rotation=0)
plt.title(r"$\langle w_p(t) w_p(t + \tau) \rangle $")
plt.xlim([time[0], time[-1]])
plt.ylim([evalZ[0], evalZ[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

imgname = imgdir + "autocorr-time-wp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

