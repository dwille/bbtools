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

# set up output file paths
infoFile = datadir + "info"
nEvenFile = datadir + "number-dens-coeffs-even"
nOddFile = datadir + "number-dens-coeffs-odd"

######################
### data analysis ####
######################

# Pull output data
# each row is a different time, transpose so column is time
nEven = np.genfromtxt(nEvenFile).T
nOdd = np.genfromtxt(nOddFile).T
nTotal = nEven + 1j*nOdd
order = np.arange(np.shape(nEven)[0])

# Calculate vfrac coefficients
vfEven = np.zeros(np.shape(nEven))
vfOdd = np.zeros(np.shape(nOdd))
vfTotal = np.zeros(np.shape(nTotal), dtype=complex)
vf_pspec = np.zeros(np.shape(nTotal))
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

# Find magnitude of coeffs and then average over time
vfMag = np.absolute(vfTotal)**2
vf_coeffs_mean = np.mean(vfMag,1)

######################
### plotting #########
######################
size=(2,2)

### mean |coeffs| over time ###
fig = plt.figure(figsize=size)
ax2 = fig.add_subplot(111)

plt.plot(order[1:], vf_coeffs_mean[1:]/vf_coeffs_mean[0], 'ko-', markerfacecolor="None",
  markersize=2.5)

ax2.set_xlabel(r"$k$")
ax2.set_ylabel(r"$\langle |\phi_k|^2 \rangle_t / \phi_0^2$")

ax2.xaxis.set_major_locator(MultipleLocator(2))
ax2.xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "vf-coeffs"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

### vfrac coeffs colormap ###
# fig = plt.figure(figsize=(3,5.5))
# ax1 = fig.add_subplot(211)
# plt.imshow(vfMag, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], order[0]-0.5, order[-1]+0.5],
#   vmin=0, vmax=0.025)
# plt.colorbar()
# plt.xlabel(r"$t\ [s]$")
# plt.ylabel(r"$order\ \ell$")
# plt.title(r'$|\phi_\ell(t)|$')
# 
# 
# # ylim for both plots
# orderMax = 20
# 
