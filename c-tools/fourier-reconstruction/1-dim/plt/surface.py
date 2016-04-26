#!/usr/bin/env python2

import sys, os
import matplotlib.pyplot as plt
import numpy as np

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
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
imgdir = root + simdir + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# set up output file paths
infoFile = datadir + "info"
nDensFile = datadir + "number-density"
vFracFile = datadir + "volume-fraction"
upFile = datadir + "part-u"
vpFile = datadir + "part-v"
wpFile = datadir + "part-w"

# Find time and evalZ
time = np.genfromtxt(infoFile, skip_footer=1)[1:]
print "      Steady state time set to: " + str(time[ts])
time = time[ts:] - time[ts]
nt = np.size(time)
evalZ = np.genfromtxt(infoFile, skip_header=1)[1:] / partR

# Find output data
numDens = np.genfromtxt(nDensFile).T
vFrac = np.genfromtxt(vFracFile).T
up = np.genfromtxt(upFile).T
vp = np.genfromtxt(vpFile).T
wp = np.genfromtxt(wpFile).T

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

## NUMBER DENSITY ##
nDensFig = plt.figure()
plt.imshow(numDens, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.colorbar(format="%.1e")
plt.title(r"$n$")
plt.xlabel(r"$t - t_0$")
plt.xticks(np.floor(np.arange(time[0], time[-1], 1000)))
plt.ylabel(r'$z/a$')

imgname = imgdir + "number-density"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## VOLUME FRACTION ##
vFracFig = plt.figure()
plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
#  vmin=0., vmax=0.5)
cbar = plt.colorbar()
plt.title(r'$Volume\ Fraction$')
plt.xlabel(r"$t - t_0$")
plt.xticks(np.floor(np.arange(time[0], time[-1], 1000)))
plt.ylabel(r'$z/a$')

imgname = imgdir + "volume-fraction"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## U_PART ##
upFig = plt.figure()
plt.imshow(up, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.colorbar()
plt.title(r"$U_p [mm/ms]$")
plt.xlabel(r"$t - t_0$")
plt.xticks(np.floor(np.arange(time[0], time[-1], 1000)))
plt.ylabel(r'$z/a$')

imgname = imgdir + "part-u"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## V_PART ##
vpFig = plt.figure()
plt.imshow(vp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.colorbar()
plt.title(r"$V_p [mm/ms]$")
plt.xlabel(r"$t - t_0$")
plt.xticks(np.floor(np.arange(time[0], time[-1], 1000)))
plt.ylabel(r'$z/a$')

imgname = imgdir + "part-v"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## W_PART ##
wpFig = plt.figure()
plt.imshow(wp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
cbar = plt.colorbar()
plt.title(r"$W_p [mm/ms]$")
plt.xlabel(r"$t - t_0$")
plt.xticks(np.floor(np.arange(time[0], time[-1], 1000)))
plt.ylabel(r'$z/a$')

imgname = imgdir + "part-w"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
