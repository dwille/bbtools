#!/usr/bin/env python2

from setup import *

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.fftpack as scifft

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

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
for tt, tval in enumerate(time):
  (vfAutoCorrReal[:,tt], vfAutoCorrImag[:,tt], vfPowerSpec[:,tt]) \
    = AutoCorrelationSpaceFFT(vFrac[:,tt])
  #vfPowerSpec[:,tt] = np.absolute(scifft.fft(vfAutoCorrReal[:,tt]))**2
 
vfPowerSpec /= (nz*nz)

# Make sure imag part is negligible
if (np.max(vfAutoCorrImag) > 5e-15) or (np.min(vfAutoCorrImag) < -5e-15):
  print "Imaginary part is not numerical zero:"
  print "   Max = %.2e" % np.max(vfAutoCorrImag)
  print "   Min = %.2e" % np.min(vfAutoCorrImag)
  sys.exit()

size = (2,2)
# # Plot # #
fig1 = plt.figure(figsize=(2,2))

# Plot MEAN of all power spectrums averaged over tt
vfAvgPowerSpec = np.mean(vfPowerSpec, 1)

ax1 = fig1.add_subplot(111)
plt.semilogx(vfAvgPowerSpec[1:]/np.max(vfAvgPowerSpec), np.arange(1,np.size(freq)), 
  'ko-', markerfacecolor="None", markersize=2.5)

ax1.set_xlabel(r'$Power/Power_{max}$')
#ax1.set_xlim([0, 1.25])
#ax1.set_xticks([0, 0.5, 1., 1.25])
ax1.set_xlim([1e-3, 1.25])
ax1.set_ylabel(r"$\ell$",rotation=0)
ax1.set_ylim([0, 15])

xy=get_axis_limits(ax1)
ax1.annotate(r"$(b)$", (0.3, 12.75))

#ax1.xaxis.set_major_locator(MultipleLocator(0.5))
#ax1.xaxis.set_minor_locator(MultipleLocator(0.25))

# SAVE
imgname = imgdir + "avg-power-spectrum-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Plot mean of autocorrelation
fig2 = plt.figure(figsize=(2,2))
ax2 = fig2.add_subplot(111)

vfAutoCorrMean = np.mean(vfAutoCorrReal, 1)
plt.plot(vfAutoCorrMean, dz, 'k-')

ax2.set_xlabel(r'$R(\Delta z)$',rotation=0)
ax2.set_xlim([-0.5, 1.])
ax2.set_ylabel(r"$\Delta z$")
ax2.set_ylim([dz[0], dz[-1]])

ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.xaxis.set_minor_locator(MultipleLocator(0.25))

imgname = imgdir + "avg-autocorr-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

#vfFig = plt.figure(figsize=(6,4))
#gs = matplotlib.gridspec.GridSpec(2,4)
#
## Plot autocorrelation (vertical slices at constant time)
#ax1 = vfFig.add_subplot(gs[0,0:2])
#plt.imshow(vfAutoCorrReal, origin="lower", aspect="auto", interpolation="none",
#  extent=[time[0], time[-1], dz[0], dz[-1]],
#  vmin=-1., vmax=1., cmap='seismic')
#
## Plot power spectrum as a colormap
#ax2 = vfFig.add_subplot(gs[0,2:4])
#plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0]-0.5*dt, time[-1]+0.5*dt, 0-0.5, np.size(freq)+0.5],
#   vmin=-np.max(vfPowerSpec), vmax=np.max(vfPowerSpec), cmap='seismic')
#
#plt.ylim([-0.5, 15.5])
#plt.ylabel(r"$k$")
# 
#imgname = imgdir + "power-spectrum-space-vf"
#plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
