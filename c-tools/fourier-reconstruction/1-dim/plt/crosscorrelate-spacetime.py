#!/usr/bin/env python2

import sys, os
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
zs = 0

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
print "      Starting z index set to: " + str(zs)

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
print "      Steady state time set to: " + str(time[ts])
time = time[ts:] - time[ts]
nt = np.size(time)

evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] # / partR
print "      Starting z set to: " + str(evalZ[zs])
nz = np.size(evalZ)
dz = evalZ - evalZ[zs]

# Find output data -- each column is a different time
numDens = np.genfromtxt(nDensFile).T[:,ts:]
vFrac = np.genfromtxt(vFracFile).T[:,ts:]
up = np.genfromtxt(upFile).T[:,ts:]
vp = np.genfromtxt(vpFile).T[:,ts:]
wp = np.genfromtxt(wpFile).T[:,ts:]

# Plot specs
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
#plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=11, numpoints=3)
plt.rc('lines', markersize=2)
plt.rc('savefig', dpi=250)
labelx = -0.17

# Crosscorrelation of volume fraction
vfFig = plt.figure()
vfCrossCorr = np.zeros((nz, nt))
vfMaxima = np.zeros((nz,3))
vfFirstMaxima = np.zeros((nz,3))
for zz, zval in enumerate(evalZ):
  if zs + zz >= nz:
    zInd = zs + zz - nz
  else:
    zInd = zs + zz
  # length of result is ceil(length(time)/2)
  vfCrossCorr[zz,:] = CrossCorrelationFFT(vFrac[zs,:], vFrac[zInd,:])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(vfCrossCorr[zz,:]))) < 0).nonzero()[0] + 1
  maxima = vfCrossCorr[zz,maximaLoc]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfMaxima[zz,0] = np.nan
    vfMaxima[zz,1] = dz[zz]
    vfMaxima[zz,2] = np.nan
  else:
    vfMaxima[zz,0] = time[maximaLoc[maxInd]]
    vfMaxima[zz,1] = dz[zz]
    vfMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[maxInd]]

    vfFirstMaxima[zz,0] = time[maximaLoc[0]]
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[0]]

vfCrossCorr /= vfCrossCorr[0,0]

## Find slope of maxima
for zz in np.arange(1,nz):
  # if the time changes drastically, make not of it
  if np.abs(vfFirstMaxima[zz,0] - vfFirstMaxima[zz-1,0]) > .1:
    firstWave = vfFirstMaxima[0:zz-1,:]

tau = firstWave[:,0]
tau = tau[:,np.newaxis]
dzMax = firstWave[:,1]
dzMax = dzMax[:,np.newaxis]

# Make sure that x,y[0] = 0
tau[0] = 0
dzMax[0] = 0

# Fit curve, assume 0 intercept
p, _, _, _ = np.linalg.lstsq(tau, dzMax)
xFit = firstWave[:,0]
yFit = p[0,0]*xFit

# Plot
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
plt.ylim([0, dz[-1]])

plt.title(r"$\langle \alpha(t,z) \alpha(t + \tau, z + \Delta z) \rangle,\ zs ="+
  str(zs) + "$")

imgname = imgdir + "crosscorr-spacetime-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Crosscorrelation of wp -- time
wpFig = plt.figure()
wpCrossCorr = np.zeros((nz, nt))
for zz, zval in enumerate(evalZ):
  if zs + zz >= nz:
    zInd = zs + zz - nz
  else:
    zInd = zs + zz
  # length of result is ceil(length(time)/2)
  wpCrossCorr[zz,:] = CrossCorrelationFFT(wp[zs,:],wp[zInd,:])

wpCrossCorr /= wpCrossCorr[0,0]


plt.imshow(wpCrossCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()

plt.xlabel(r"$\tau\ [s]$")
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
plt.ylabel(r"$dz\ [mm]$", rotation=0)
plt.title(r"$\langle w_p(t,z) w_p(t + \tau, z + \Delta z) \rangle,\ zs ="+
  str(zs) + "$")

imgname = imgdir + "crosscorr-spacetime-wp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
