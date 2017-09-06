#!/usr/bin/env python2
from setup import *
os.system("clear")

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Power Spectrum -- Time"
print ""

######################
### initialization ###
######################

# Get simulation parameters
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

# Find output data -- each column is a different time step
vFracFile = datadir + "volume-fraction"
vFracPull = np.genfromtxt(vFracFile).T[:,tsInd:]

# Interpolate data to constant-dt time
interp_dt = np.mean(np.diff(time))
uniform_time = np.arange(0, interp_dt*np.size(time), interp_dt)
nt = np.size(uniform_time)

vFrac = np.zeros((nz, nt))
for zz in np.arange(0, nz):
  vFrac[zz,:] = np.interp(uniform_time, time, vFracPull[zz,:])

time = uniform_time
dt = interp_dt

# Fourier Transform Parameters
freq = scifft.fftfreq(nt, dt)

# Autocorrelation of volume fraction, power spectrum of autocorrelation
vfAutoCorr = np.zeros((nz, nt))
vfPowerSpec = np.zeros((nz,nt))
for zz, zval in enumerate(evalZ):
  # length of result is ceil(length(time)/2)
  vfAutoCorr[zz,:] = AutoCorrelationFFT(vFrac[zz,:])
  vfPowerSpec[zz,:] = np.absolute(scifft.fft(vfAutoCorr[zz,:]))**2

vfPowerSpec /= (nt*nt)

# Find mean of autocorr, pspec (avg over z)
vfAutoCorrMean = np.mean(vfAutoCorr, 0)
vfPowerSpecMean = np.mean(vfPowerSpec, 0)

# Find most prominant frequency based on max power
freqMax = freq[np.argmax(vfPowerSpecMean)]
powerMax = np.max(vfPowerSpecMean)

######################
### plotting #########
######################
size = (2,2)

### Power Spectrum ###
fig1 = plt.figure(figsize=size)
ax1 = fig1.add_subplot(111)

plt.plot(freq/freqMax, vfPowerSpecMean/powerMax, 'ko-', markerfacecolor="None",
  markersize=2.5)

ax1.set_xlabel(r'$f/f_{max}$')
ax1.set_xlim([0, 7])

ax1.set_ylabel(r'$P/P_{max}$')
ax1.set_ylim([0, 1.25])
ax1.set_yticks([0,0.5,1,1.25])
ax1.yaxis.set_minor_locator(MultipleLocator(0.25))

#firstThreeMaxIndices = np.argsort(-vfPowerSpecMean[0:len(freq)/2])[0:3]
plt.text(freqMax + 0.4, 1 + 0.05, r"$f_{max} = %.4f\ [\mathrm{Hz}]$" % freqMax)

# save
imgname = imgdir + "avg-power-spectrum-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')


### Autocorrelation ###
fig2 = plt.figure(figsize=size)
ax2 = fig2.add_subplot(111)

plt.plot(time, vfAutoCorrMean, 'k-')

ax2.set_xlim([0., time[-1]])
ax2.set_xlabel(r"$\tau\ [\mathrm{s}]$")
ax2.set_ylim([-0.5, 1])
ax2.set_ylabel(r"$R(\tau)$")
for i in np.arange(1,5):
  xpoints = [i/freqMax,i/freqMax]
  ypoints = [-1,1]
  plt.plot(xpoints, ypoints, 'k:')

# save
imgname = imgdir + "avg-autocorr-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "\n      ...Done!"

# colormap plots
#  # Plot
#  vfFig = plt.figure(figsize=(3,2))
#  #gs = matplotlib.gridspec.GridSpec(2,4)
#  
#  # Plot auto-correlation (slices at constant location)
#  ax1 = vfFig.add_subplot(gs[0,0:2])
#  plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
#    extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
#    vmin=-1,vmax=1,cmap='seismic')
#  
#  ax1.set_xlabel(r'$Time\ [s]$')
#  ax1.set_xlim([0, time[-1]])
#  ax1.set_ylabel(r'$z\ [mm]$')
#  ax1.set_ylim([np.min(evalZ), np.max(evalZ)])
#  
#  # Plot ALL power spectrums as colormap
#  ax2 = vfFig.add_subplot(gs[0,2:4])
#  plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
#    extent=[0, np.max(freq*2), evalZ[0], evalZ[-1]],
#    vmin=-np.max(vfPowerSpec), vmax=np.max(vfPowerSpec), cmap='seismic')
#  
#  ax2.set_xlim([0, 7])
#  ax2.xaxis.set_ticklabels([])
#  ax2.set_ylim([np.min(evalZ), np.max(evalZ)])
#  ax2.yaxis.set_ticklabels([])
