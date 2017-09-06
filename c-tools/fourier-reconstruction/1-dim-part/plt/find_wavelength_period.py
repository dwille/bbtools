#!/usr/bin/env python2
from setup import *
os.system('clear')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# Setup simulation parameters
(partR, nparts, rho, vFracMean, simdir, tstart) = simParams(sys)
nu = .01715   ## mm^2/ms XXX

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# get avg volume fraction
domX = 42 # XXX
domY = 42 # XXX
domZ = 126 # XXX
avgVolumeFraction = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)

# non-dimensionalize axes (time, position)
evalZ /= (2.*partR)               ## mm / mm
tau = (2.*partR)*(2.*partR)/nu    ## mm^2/(mm^2/ms) = ms
tau /= 1000.                      ## ms -> s
#time /= tau                       ## s / s

# Interpolate time to constant dt
interp_dt = np.mean(np.diff(time))
uniform_time = np.arange(0, interp_dt*np.size(time), interp_dt)
nt = np.size(uniform_time)

# Pull Volume fraction
vfrac_pull = np.genfromtxt(datadir + "volume-fraction").T[:,tsInd:]

# Interpolate volume fraction to constant dt time
vfrac = np.zeros((nz, nt))
for zz in np.arange(0, nz):
  vfrac[zz,:] = np.interp(uniform_time, time, vfrac_pull[zz,:])

## PERIOD/FREQUENCY ##
# Loop over slicies at constant z -> phi(t)
vf_fft = np.zeros((nz, nt))
for zz in np.arange(0, nz):
  flucts = vfrac[zz,:] - np.mean(vfrac[zz,:])
  vf_fft[zz,:] = scifft.fft(flucts)

# Normalize -- calc power spectrum (if n = 2) XXX
vf_fft /= nt*nt
vf_fft = np.mean(vf_fft**2, 0)

# Get frequencies
freq = scifft.fftfreq(nt, interp_dt)

# Get max freq, power spectrum
freqMax = freq[np.argmax(vf_fft)]
powerMax = np.max(vf_fft)

# Plot
fig = plt.figure(figsize = (2,2))
ax = fig.add_subplot(111)

plt.plot(freq/freqMax, vf_fft/powerMax, 'ko-', markerfacecolor="None", markersize=2.5)
ax.set_xlabel(r'$f/f_{max}$')
ax.set_xlim([0, 7])
ax.set_ylabel(r'$P/P_{max}$')
ax.set_ylim([0, 1.25])
ax.set_yticks([0,0.5,1,1.25])
ax.yaxis.set_minor_locator(MultipleLocator(0.25))
ax.text(freqMax + 0.4, 1 + 0.05, r"$f_{max} = %.4f\ [\mathrm{Hz}]$" % freqMax)


plt.savefig(imgdir + "fft-time-ps.png", bbox_inches='tight', format='png')

## WAVELENGTH ##
# Loop over slices at constant t -> phi(z)
#for tt,tval in enumerate(time):

  

   

