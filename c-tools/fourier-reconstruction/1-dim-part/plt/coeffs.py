#!/usr/bin/env python2

import sys, os
import matplotlib.pyplot as plt
import numpy as np

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                  Coefficients"
print ""

# SIMULATION PARAMETERS
partR = 2.1
domX = 42
domY = 42
domZ = 126
nparts = 2000

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
nEvenFile = datadir + "number-dens-coeffs-even"
nOddFile = datadir + "number-dens-coeffs-odd"
wpEvenFile = datadir + "part-w-coeffs-even"
wpOddFile = datadir + "part-w-coeffs-odd"

# Find time and evalZ
time = np.genfromtxt(infoFile, skip_footer=1)[1:]
evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] / partR

# Find output data -- each column is a timestep
nEven = np.genfromtxt(nEvenFile).T
nOdd = np.genfromtxt(nOddFile).T
wpEven = np.genfromtxt(wpEvenFile).T
wpOdd = np.genfromtxt(wpOddFile).T
order = np.arange(np.shape(wpEven)[0])

# Calculate vfrac
vfEven = np.zeros(np.shape(nEven))
vfOdd = np.zeros(np.shape(nOdd))
for oo in order:
  if oo == 0:
    base = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)
    vfEven[oo,:] = 0.5*base
    vfOdd[oo,:] = 0.5*base
  elif oo != 0:
    k = 2.*np.pi*oo/domZ
    ka = k*partR
    correction = 4.*np.pi/(k*k*k)*(np.sin(ka) - ka*np.cos(ka))
    vfEven[oo,:] = correction*nEven[oo,:]
    vfOdd[oo,:] = correction*nOdd[oo,:]

# Find magnitude of coeffs
wpCoeffs = 0.5*np.sqrt(np.square(wpEven) + np.square(wpOdd));
vfCoeffs = 0.5*np.sqrt(np.square(vfEven) + np.square(vfOdd));

# Plot specs
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
#plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=11, numpoints=3)
plt.rc('lines', markersize=4)
plt.rc('savefig', dpi=250)
labelx = -0.17

## vfrac coeffs ##
vFracFig = plt.figure()
plt.imshow(vfCoeffs, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], order[0], order[-1]],
  vmin=0, vmax=0.01)
plt.colorbar()
plt.xlabel('time')
plt.ylabel('order')
plt.title(r'$|\beta_\ell|$')

imgname = imgdir + "vf-coeffs"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## wp-coeffs ##
wpFig = plt.figure()
plt.imshow(wpCoeffs, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], order[0], order[-1]])
plt.colorbar()
plt.xlabel('time')
plt.ylabel('order')
plt.title(r'$|(nw)_\ell|$')

imgname = imgdir + "wp-coeffs"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
