#!/usr/bin/env python
#
# Calculate the spatial power spectrum of the volume fraction.
# 1) Take the reconstructed volume fraction at time t_i, phi(z, t_i). If you'd
#     like, this is a vertical "slice" in the z-t volume fraction plot.
# 2) Subtract the mean of phi(z, t_i), which is just the bulk volume fraction.
#     Let's call it phi_prime = phi - phi_mean
# 3) Take the DFT of the resulting signal phi_prime to get phi_hat,
#     phi_hat = DFT(phi_prime)
# 4) Square phi_hat to get the power spectrum
# 5) Repeat steps 1-4 for all t_i, then average over t.
# 
# Some notes:
# a) I'm no longer sure why I subtracted the mean in step 2. I'm doing this
#     calculation at the same time as I perform an autocorrelation, which
#     traditionally uses zero-mean signals. I don't think this affects the power
#     spectrum outside of the zero-wavenumber term since the FT of phi_mean
#     (which is constant) is a delta function.
# 
# b) The square in step 4 occurs after the DFT rather than before. This is
#     termed the "Energy spectral density" from [1] (see equation for S_xx);
#     I've also seen this called the "power spectrum" which is the term I've
#     been using.
# 
# c) The "l" on the x-axis represents the discrete frequencies (wavenumber,
#     technically, since we're looking at space) from the DFT, or the "cycles
#     per Nz samples" (Nz being the number of grid cells in the vertical
#     direction, or the length of the discrete signal phi(z, t_i) ). This does
#     NOT use the convention with 2pi. The wavelength is lambda = (Nz * dz)/l or
#     lambda = Lz/l. So, l = 0 is a constant wave (infinite wavelength, the avg
#     volume fraction), l = 1 corresponds to a wave with a wavelength = Lz,
#     l = 2 corresponds to a wave with wavelength = Lz / 2
#
#  We can find the wavelengths by using this:
#  print(1./np.fft.fftfreq(nz, np.mean(np.diff(evalZ))))
# 

from setup import *

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.fftpack as scifft

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

os.system("clear")

######################
### initialization ###
######################

# Setup simulation parameters
(partR, nparts, rho, vFracMean, simdir, tstart) = simParams(sys)
nu = .01715   ## mm^2/ms XXX

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

######################
### data analysis ####
######################

# Find output data -- each row is a different time, tranpose so column is time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Autocorrelation of volume fraction (space), power spectrum of result
vfAutoCorr = np.zeros((nz,nt))
vfPowerSpec = np.zeros((nz,nt))
for tt, tval in enumerate(time):
  (vfAutoCorr[:,tt], vfPowerSpec[:,tt]) = AutoCorrelationSpaceFFT(vFrac[:,tt])

  #vfAutoCorr[:,tt] = AutoCorrelationFFT(vFrac[:,tt])
  #vfPowerSpec[:,tt] = np.absolute(scifft.fft(vfAutoCorr[:,tt]))**2
 
vfPowerSpec /= (nz*nz)
# XXX how do the wavelength lie with the wrap-around frequencies??

# Find mean of autocorr, pspec (avg over t)
vfAutoCorrMean = np.mean(vfAutoCorr, 1)
vfPowerSpecMean = np.mean(vfPowerSpec, 1)

# Fourier transform parameters
dz = np.mean(np.diff(evalZ - evalZ[0]))
inv_wlength = scifft.fftfreq(nz, dz) # 1/wavelength
modes = np.arange(0, np.size(inv_wlength))

# Find most prominant wavelength/mode based on max power
max_mode = modes[np.argmax(vfPowerSpecMean)]
max_power = np.max(vfPowerSpecMean)

######################
### plotting #########
######################
size = (2,2)

### Power Spectrum ###
fig1 = plt.figure(figsize=(2,2))
ax1 = fig1.add_subplot(111)
plt.plot(modes, vfPowerSpecMean/max_power, 'ko-', markerfacecolor="None",
  markersize=2.5)

ax1.set_xlabel(r"$L_z / \lambda$", fontsize=14)
ax1.set_xlim([0, 15])

ax1.set_ylabel(r'$P/P_{max}$', fontsize=14)
ax1.set_ylim([0, 1.25])

ax1.xaxis.set_major_locator(MultipleLocator(5))
ax1.xaxis.set_minor_locator(MultipleLocator(1))
ax1.yaxis.set_major_locator(MultipleLocator(0.50))
ax1.yaxis.set_minor_locator(MultipleLocator(0.25))

#ax1.set_yticks([0,0.5,1,1.25])

for tick in ax1.xaxis.get_major_ticks():
  tick.label.set_fontsize(12)
for tick in ax1.yaxis.get_major_ticks():
  tick.label.set_fontsize(12)

# SAVE
imgname = imgdir + "avg-power-spectrum-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')

# Plot mean of autocorrelation
fig2 = plt.figure(figsize=(2,2))
ax2 = fig2.add_subplot(111)

plt.plot(evalZ, vfAutoCorrMean, 'k-')

ax2.set_xlabel(r"$z\ [mm]$")
ax2.set_xlim([evalZ[0], evalZ[-1]])

ax2.set_ylabel(r'$R(z)$',rotation=0)
ax2.set_ylim([-0.5, 1.])

ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(MultipleLocator(0.25))

imgname = imgdir + "avg-autocorr-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

#  #vfFig = plt.figure(figsize=(6,4))
#  #gs = matplotlib.gridspec.GridSpec(2,4)
#  #
#  ## Plot autocorrelation (vertical slices at constant time)
#  #ax1 = vfFig.add_subplot(gs[0,0:2])
#  #plt.imshow(vfAutoCorrReal, origin="lower", aspect="auto", interpolation="none",
#  #  extent=[time[0], time[-1], dz[0], dz[-1]],
#  #  vmin=-1., vmax=1., cmap='seismic')
#  #
#  ## Plot power spectrum as a colormap
#  #ax2 = vfFig.add_subplot(gs[0,2:4])
#  #plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
#  #   extent=[time[0]-0.5*dt, time[-1]+0.5*dt, 0-0.5, np.size(wlength)+0.5],
#  #   vmin=-np.max(vfPowerSpec), vmax=np.max(vfPowerSpec), cmap='seismic')
#  #
#  #plt.ylim([-0.5, 15.5])
#  #plt.ylabel(r"$k$")
#  # 
#  #imgname = imgdir + "power-spectrum-space-vf"
#  #plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Make sure imag part is negligible
##if (np.max(vfAutoCorrImag) > 5e-15) or (np.min(vfAutoCorrImag) < -5e-15):
##  print "Imaginary part is not numerical zero:"
##  print "   Max = %.2e" % np.max(vfAutoCorrImag)
##  print "   Min = %.2e" % np.min(vfAutoCorrImag)
##  sys.exit()

