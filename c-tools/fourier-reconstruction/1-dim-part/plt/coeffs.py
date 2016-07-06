#!/usr/bin/env python2
from setup import *
os.system('clear')

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                  Coefficients"
print ""

# SIMULATION PARAMETERS
domX = 42
domY = 42
domZ = 126
nparts = 2000

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# set up output file paths
infoFile = datadir + "info"
nEvenFile = datadir + "number-dens-coeffs-even"
nOddFile = datadir + "number-dens-coeffs-odd"

# Find output data -- each column is a timestep
nEven = np.genfromtxt(nEvenFile).T
nOdd = np.genfromtxt(nOddFile).T
nTotal = nEven + 1j*nOdd
order = np.arange(np.shape(nEven)[0])

# Calculate vfrac
vfEven = np.zeros(np.shape(nEven))
vfOdd = np.zeros(np.shape(nOdd))
vfTotal = np.zeros(np.shape(nTotal), dtype=complex)
for oo in order:
  if oo == 0:
    base = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)
    vfEven[oo,:] = 0.5*base
    vfOdd[oo,:] = 0.5*base
    vfTotal[oo,:] = base
  elif oo != 0:
    k = 2.*np.pi*oo/domZ
    ka = k*partR
    correction = 4.*np.pi/(k*k*k)*(np.sin(ka) - ka*np.cos(ka))
    vfEven[oo,:] = correction*nEven[oo,:]
    vfOdd[oo,:] = correction*nOdd[oo,:]
    vfTotal[oo,:] = correction*nTotal[oo,:]

# Find magnitude of coeffs
vfMag = np.absolute(vfTotal)

fig = plt.figure()
## vfrac coeffs ##
ax1 = fig.add_subplot(211)
plt.imshow(vfMag, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], order[0]-0.5, order[-1]+0.5],
  vmin=0, vmax=0.025)
plt.colorbar()
plt.xlabel('time')
plt.ylabel('order')
plt.title(r'$|\phi_\ell|$')

# mean
vf_coeffs_mean = np.mean(vfMag,1)
vf_coeffs_sdev = np.std(vfMag,1)
ax2 = fig.add_subplot(212)
plt.semilogx(vf_coeffs_mean, order, 'ko--')
#plt.semilogx(vf_coeffs_mean + vf_coeffs_sdev, order, 'ro:')
#plt.semilogx(vf_coeffs_mean - vf_coeffs_sdev, order, 'ro:')

ax2.set_xlim([5e-4, 5e-1])

#ax2.xaxis.set_major_locator(MultipleLocator(.25))
#ax2.xaxis.set_minor_locator(MultipleLocator(.05))
ax2.yaxis.set_major_locator(MultipleLocator(5))
ax2.yaxis.set_minor_locator(MultipleLocator(1))
ax2.grid(True)

imgname = imgdir + "vf-coeffs"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
