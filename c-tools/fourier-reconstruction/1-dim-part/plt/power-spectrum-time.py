#!/usr/bin/env python2
from setup import *
os.system("clear")

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Power Spectrum -- Time"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
#(root, simdir, datadir, imgdir) = directoryStructureDevel(simdir)
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFracPull = np.genfromtxt(vFracFile).T[:,tsInd:]

# Interpolate data to constant dt time
interpDT = np.mean(np.diff(time))
uniform_time = np.arange(0, interpDT*np.size(time), interpDT)
nt = np.size(uniform_time)

vFrac = np.zeros((nz, nt))
for zz in np.arange(0, nz):
  vFrac[zz,:] = np.interp(uniform_time, time, vFracPull[zz,:])

time = uniform_time

# Fourier Transform Parameters
dt = interpDT
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

size = (2,2)
# PLOT #
fig1 = plt.figure(figsize=size)

# Plot MEAN of all power spectrums averaged over zz
vfAvgPowerSpec = np.mean(vfPowerSpec, 0)
freqMax = freq[np.argmax(vfAvgPowerSpec)]
powerMax = np.max(vfAvgPowerSpec)

ax3 = fig1.add_subplot(111)
plt.plot(freq/freqMax, vfAvgPowerSpec/powerMax, 'ko-', markerfacecolor="None", markersize=2.5)

ax3.set_xlabel(r'$f/f_{max}\ [\mathrm{Hz}]$')
ax3.set_xlim([0, 7])
ax3.set_ylabel(r'$Power/P_{max}$')
ax3.set_ylim([0, 1.25])
ax3.set_yticks([0,0.5,1,1.25])

ax3.yaxis.set_minor_locator(MultipleLocator(0.25))

ax3.annotate(r"$(d)$",xy=get_axis_limits(ax3))

firstThreeMaxIndices = np.argsort(-vfAvgPowerSpec[0:len(freq)/2])[0:3]
powerMax = vfAvgPowerSpec[np.argmax(vfAvgPowerSpec)]
plt.text(freqMax + 0.4, 1, r"$f_{max} = %.2f$" % freqMax)

# SAVE
imgname = imgdir + "avg-power-spectrum-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

fig2 = plt.figure(figsize=size)
# Plot mean of autocorrelation
ax4 = fig2.add_subplot(111)
plt.plot(time, vfAutoCorrMean, 'k-')

ax4.set_xlim([0., time[-1]])
ax4.set_xlabel(r"$\tau\ [\mathrm{s}]$")
ax4.set_ylim([-0.5, 1])
ax4.set_ylabel(r"$R(\tau)$")
for i in np.arange(1,5):
  xpoints = [i/freqMax,i/freqMax]
  ypoints = [-1,1]
  plt.plot(xpoints, ypoints, 'k:')

ax4.annotate(r"$(b)$",xy=get_axis_limits(ax4))

# SAVE
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
