#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "             Avg Crosscorrelation Plot"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
#(root, simdir, datadir, imgdir) = directoryStructureDevel(simdir)
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Pull output data -- each column is a different time
vfSaveFile = datadir + "z-averaged-vfrac-xcorr"
vfCrossCorr = np.genfromtxt(vfSaveFile)

# Find maxima
vfFirstMaxima = np.zeros((nz,3))
vfSecondMaxima = np.zeros((nz,3))
vfThirdMaxima = np.zeros((nz,3))
for zz in np.arange(0,nz):
  # Find maxima locations by using where second derivative changes
  maximaLoc = (np.diff(np.sign(np.diff(vfCrossCorr[zz,:]))) < 0).nonzero()[0] + 1
  # maximum values
  maxima = vfCrossCorr[zz,maximaLoc]
  # maximum indices
  maxInd = np.argmax(maxima)

  if np.size(maximaLoc) == 0:
    vfFirstMaxima[zz,0] = np.nan
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = np.nan

    vfSecondMaxima[zz,0] = np.nan
    vfSecondMaxima[zz,1] = dz[zz]
    vfSecondMaxima[zz,2] = np.nan
  elif np.size(maximaLoc) > 0:
    vfFirstMaxima[zz,0] = time[maximaLoc[0]]
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[0]]

    if np.size(maximaLoc) > 1:
      vfSecondMaxima[zz,0] = time[maximaLoc[1]]
      vfSecondMaxima[zz,1] = dz[zz]
      vfSecondMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[1]]

      if np.size(maximaLoc) > 2:
        vfThirdMaxima[zz,0] = time[maximaLoc[2]]
        vfThirdMaxima[zz,1] = dz[zz]
        vfThirdMaxima[zz,2] = vfCrossCorr[zz,maximaLoc[2]]

  if vfCrossCorr[zz,0] > vfCrossCorr[zz,1]:
    if np.size(maximaLoc) > 2:
      vfThirdMaxima[zz,0] = vfSecondMaxima[zz,0]
      vfThirdMaxima[zz,1] = vfSecondMaxima[zz,1]
      vfThirdMaxima[zz,2] = vfSecondMaxima[zz,2]

    if np.size(maximaLoc) > 1:
      vfSecondMaxima[zz,0] = vfFirstMaxima[zz,0]
      vfSecondMaxima[zz,1] = vfFirstMaxima[zz,1]
      vfSecondMaxima[zz,2] = vfFirstMaxima[zz,2]

    vfFirstMaxima[zz,0] = time[0]
    vfFirstMaxima[zz,1] = dz[zz]
    vfFirstMaxima[zz,2] = vfCrossCorr[zz,0]

# Find slope of maxima -- wavespeed
tau = -np.ones(nz)
zeta = -np.ones(nz)
tau[0] = vfFirstMaxima[0,0]
zeta[0] = vfFirstMaxima[0,1]
for zz in np.arange(0,nz-1):
  tauFirstNext = vfFirstMaxima[zz+1,0]
  tauSecondNext = vfSecondMaxima[zz+1,0]
  tauThirdNext = vfThirdMaxima[zz+1, 0]

  # Use first max if it is close
  if np.abs(tau[zz] - tauFirstNext) < 0.01:
    tau[zz+1] = tauFirstNext
    zeta[zz+1] = vfFirstMaxima[zz+1,1]
  # if it's too far, try the second max
  elif np.abs(tau[zz] - tauSecondNext) < 0.01:
    tau[zz+1] = tauSecondNext
    zeta[zz+1] = vfSecondMaxima[zz+1,1]
  # if it's too far, try the third
  elif np.abs(tau[zz] - tauThirdNext) < 0.01:
    tau[zz+1] = tauThirdNext
    zeta[zz+1] = vfThirdMaxima[zz+1,1]
  # if that is not good, quit
  else:
    break;

tau = tau[tau != -1]
zeta = zeta[zeta != -1]
  
# Fit curve, assume 0 intercept -- p is slope
p, _, _, _ = np.linalg.lstsq(tau[:,np.newaxis], zeta)
yFit = p*tau

## PLOTTING ##
# Volume Fraction
vfFig = plt.figure()
plt.imshow(vfCrossCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()

#plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'wo', alpha=0.4)
#plt.plot(vfSecondMaxima[:,0], vfSecondMaxima[:,1], 'ro', alpha=0.4)
#plt.plot(vfThirdMaxima[:,0], vfThirdMaxima[:,1], 'go', alpha=0.4)
plt.plot(tau, zeta, 'k-')

plt.plot(tau, yFit, 'w--')
cTxtString = r"$dz = %.4f\tau$" % p
#plt.text(time[len(time)/2], dz[len(dz)/2], cTxtString, fontsize=12)

plt.xlabel(r"$\tau\ [s]$")
plt.xlim([0, time[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

plt.ylabel(r"$\Delta z\ [mm]$", rotation=0, labelpad=20)
plt.ylim([dz[0], dz[-1]])

plt.title(r"$\langle \phi(t,z) \phi(t + \tau, z + \Delta z) \rangle$")

imgname = imgdir + "avg-crosscorr-spacetime-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# # Particle vertical velocity
# wpFig = plt.figure()
# plt.imshow(wpCrossCorr, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], dz[0], dz[-1]])
# plt.colorbar()
# 
# plt.xlabel(r"$\tau\ [s]$")
# plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
# plt.ylabel(r"$dz\ [mm]$", rotation=0)
# plt.title(r"$\langle w_p(t,z) w_p(t + \tau, z + \Delta z) \rangle$")
# 
# imgname = imgdir + "avg-crosscorr-spacetime-wp"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
print "      Done!"
