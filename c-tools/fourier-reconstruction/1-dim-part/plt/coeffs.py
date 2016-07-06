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

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

nparts = float(simdir.partition('/')[0])

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

fig = plt.figure(figsize=(3,5.5))
## vfrac coeffs ##
ax1 = fig.add_subplot(211)
plt.imshow(vfMag, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], order[0]-0.5, order[-1]+0.5],
  vmin=0, vmax=0.025)
plt.colorbar()
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r"$order\ \ell$")
plt.title(r'$|\phi_\ell(t)|$')

# cesaro mean
for oo in order:
  vfMag[oo] *= (1. - oo/(order[-1] + 1.))
# mean
vf_coeffs_mean = np.mean(vfMag,1)
vf_coeffs_med = np.median(vfMag,1)
vf_coeffs_max = np.max(vfMag,1)
vf_coeffs_sdev = np.std(vfMag,1)

ax2 = fig.add_subplot(212)
plt.semilogx(vf_coeffs_mean, order, 'ko:')
plt.semilogx(vf_coeffs_med, order, 'ro:')
plt.semilogx(vf_coeffs_mean + vf_coeffs_sdev, order, 'k--')
plt.semilogx(vf_coeffs_mean - vf_coeffs_sdev, order, 'k--')
plt.semilogx(vf_coeffs_max, order, 'bo:')

ax2.set_ylabel(r"$order\ \ell$")
ax2.set_xlabel(r"$\langle \left(1 - \frac{\ell}{L+1}\right) |\phi_\ell| \rangle_t$")
#ax2.set_xlim([.1, 20])

#ax2.xaxis.set_major_locator(MultipleLocator(.25))
#ax2.xaxis.set_minor_locator(MultipleLocator(.05))
ax2.yaxis.set_major_locator(MultipleLocator(5))
ax2.yaxis.set_minor_locator(MultipleLocator(1))
ax2.grid(True)

## std
#ax3 = fig.add_subplot(313)
#plt.semilogx(vf_coeffs_sdev/vf_coeffs_mean, order, 'ro:')
#ax3.set_xlim([.1,1])
#ax3.grid(True)



imgname = imgdir + "vf-coeffs"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
