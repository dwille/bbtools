#!/usr/bin/env python2
from setup import *
os.system("clear")

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
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Interpolate data to constant dt time
interpDT = np.mean(np.diff(time))
uniform_time = np.arange(0, interpDT*np.size(time), interpDT)
for zz in np.arange(0, nz):
  vFrac[zz,:] = np.interp(uniform_time, time, vFrac[zz,:])

time = uniform_time
dt = interpDT

# Fourier Transform Parameters
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

# Plot
plt.rc('figure', figsize=(4,6))
vfFig = plt.figure()

# Plot auto-correlation
vfFig.add_subplot(311)
plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.xlabel(r'$Time [s]$')

# Plot ALL power spectrums as colormap
vfFig.add_subplot(312)
plt.imshow(vfPowerSpec, origin="lower", aspect="auto", interpolation="none",
  extent=[0, np.max(freq*2), evalZ[0], evalZ[-1]])
plt.xlim([0, 7])

# Plot MEAN of all power spectrums averaged over zz
vfAvgPowerSpec = np.mean(vfPowerSpec, 0)

vfFig.add_subplot(313)
plt.plot(freq, vfAvgPowerSpec, 'ko:')
plt.xlabel(r'$f\ [Hz]$')
plt.ylabel(r'$PS$')
plt.xlim([0, 7])

freqMax = freq[np.argmax(vfAvgPowerSpec)]
powerMax = vfAvgPowerSpec[np.argmax(vfAvgPowerSpec)]
strText = "f = %.3f" % freqMax
plt.text(freqMax + 0.1, powerMax, strText)

# SAVE
imgname = imgdir + "power-spectrum-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
