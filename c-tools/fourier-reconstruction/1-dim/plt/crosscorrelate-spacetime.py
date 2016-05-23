#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                 Crosscorrelate"
print ""

# SIMULATION PARAMETERS
zstart = 0

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
#(root, simdir, datadir, imgdir) = directoryStructureDevel(simdir)
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
zsInd = np.argwhere(evalZ >= tstart)[0]
print "      Starting z set to: %.3f [mm]" % evalZ[zsInd]
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Crosscorrelation of volume fraction
vfFig = plt.figure()
vfCrossCorr = np.zeros((nz, nt))
vfMaxima = np.zeros((nz,3))
vfFirstMaxima = np.zeros((nz,3))

print nz, nt

for zz, zval in enumerate(evalZ):
  if zsInd + zz >= nz:
    zInd = zsInd + zz - nz
  else:
    zInd = zsInd + zz
  # length of result is ceil(length(time)/2)
  vfCrossCorr[zz,:] = CrossCorrelationFFT(vFrac[zsInd,:], vFrac[zInd,:])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(vfCrossCorr[zz,:]))) < 0).nonzero()[0] + 1
  maxima = vfCrossCorr[zz,maximaLoc]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfMaxima[zz,0] = np.nan
    vfMaxima[zz,1] = dz[zz]
    vfMaxima[zz,2] = np.nan
  else:
    vfMaxima[zz,0] = time[maximaLoc[maxInd]]
    vfMaxima[zz,1] = dz[zz]
    vfMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[maxInd]]

    vfFirstMaxima[zz,0] = time[maximaLoc[0]]
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[0]]

vfCrossCorr /= vfCrossCorr[0,0]

## Find slope of maxima
# Make sure that the first maxima is 0,0
#vfFirstMaxima[0,0] = 0
# check SO #3689710 -- peak detection in a 2d array
# check SO #3686345 -- how to find the local minima of a smooth...
print vfFirstMaxima[0:10,0]
print vfMaxima[0:10,0]
#for zz in np.arange(1,nz):
#  # if the time changes drastically, make note of it 
#  if np.abs(vfFirstMaxima[zz,0] - vfFirstMaxima[zz-1,0]) > .5:
#    firstWave = vfFirstMaxima[0:zz-1,:]
#    break
#
#tau = firstWave[:,0]
#tau = tau[:,np.newaxis]
#dzMax = firstWave[:,1]
#dzMax = dzMax[:,np.newaxis]
#
## Make sure that x,y[0] = 0
#tau[0] = 0
#dzMax[0] = 0
#
## Fit curve, assume 0 intercept
#p, _, _, _ = np.linalg.lstsq(tau, dzMax)
#xFit = firstWave[:,0]
#yFit = p[0,0]*xFit

# Plot
plt.imshow(vfCrossCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()

plt.plot(vfMaxima[:,0], vfMaxima[:,1], '--', color='0.75')
plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'k--')

#plt.plot(firstWave[:,0], firstWave[:,1], 'ko', markevery=25, 
#  markerfacecolor="None")
#plt.plot(xFit, yFit, '--')
#
#cTxtString = r"$dz = %.4f\tau$" % p[0,0]
#plt.text(xFit[-1], yFit[-1], cTxtString, fontsize=12)

plt.xlabel(r"$\tau\ [s]$")
plt.xlim([0, time[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
plt.ylabel(r"$dz\ [mm]$", rotation=0)
plt.ylim([0, dz[-1]])

plt.title(r"$\langle \alpha(t,z) \alpha(t + \tau, z + \Delta z) \rangle,\ " +
  "zs = %.3f\ [mm]$" % evalZ[zsInd])

imgname = imgdir + "crosscorr-spacetime-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

## Crosscorrelation of wp -- time
#wpFig = plt.figure()
#wpCrossCorr = np.zeros((nz, nt))
#for zz, zval in enumerate(evalZ):
#  if zs + zz >= nz:
#    zInd = zs + zz - nz
#  else:
#    zInd = zs + zz
#  # length of result is ceil(length(time)/2)
#  wpCrossCorr[zz,:] = CrossCorrelationFFT(wp[zs,:],wp[zInd,:])
#
#wpCrossCorr /= wpCrossCorr[0,0]
#
#
#plt.imshow(wpCrossCorr, origin="lower", aspect="auto", interpolation="none",
#  extent=[time[0], time[-1], dz[0], dz[-1]])
#plt.colorbar()
#
#plt.xlabel(r"$\tau\ [s]$")
#plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
#plt.ylabel(r"$dz\ [mm]$", rotation=0)
#plt.title(r"$\langle w_p(t,z) w_p(t + \tau, z + \Delta z) \rangle,\ zs ="+
#  str(zs) + "$")
#
#imgname = imgdir + "crosscorr-spacetime-wp"
#plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "      Done!"
