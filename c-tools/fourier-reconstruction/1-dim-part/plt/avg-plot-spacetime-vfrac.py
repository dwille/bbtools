#!/usr/bin/env python2
from setup import *
os.system('clear')

from matplotlib.ticker import MultipleLocator
from scipy.signal import argrelextrema

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "             Avg Crosscorrelation Plot"
print ""

# Setup simulation parameters
(partR, nparts, rho, avgVolumeFraction, simdir, tstart) = simParams(sys)
nu = .01715   ## mm^2/ms XXX

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Pull output data -- each column is a different time
vfSaveFile = datadir + "z-averaged-vfrac-xcorr"
vfCrossCorr = np.genfromtxt(vfSaveFile)

# get phase averaged data
phaseVelDir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVelDir, skip_header=1)

for nn in np.arange(0,np.size(phaseVelData, 0)):
  currN = phaseVelData[nn,0]
  currRho = phaseVelData[nn,1]
  if (nparts == currN) & (rho == currRho):
    phaseVel = phaseVelData[nn,2]
phaseVel *= 1000  # mm / ms -> mm / s

#tstar = (2.*partR) / phaseVel       ## mm / (mm/s) = s
# Non-dimensionalize
evalZ /= (2.*partR)                 ## mm / mm
dz /= (2.*partR)                    ## mm / mm
tstar = (2.*partR)*(2.*partR)/nu    ## mm^2/(mm^2/ms) = ms
tstar /= 1000.                      ## ms -> s
time /= tstar                       ## s / s
xEnd = time[-1]
xEnd = 100.*np.floor(xEnd / 100.)

# Find maxima
vfFirstMaxima = np.zeros((nz,3))
vfSecondMaxima = np.zeros((nz,3))
vfThirdMaxima = np.zeros((nz,3))
for zz in np.arange(0,nz):
  # Find maxima locations by using where second derivative changes
  maximaLoc = (np.diff(np.sign(np.diff(vfCrossCorr[zz,:]))) < 0).nonzero()[0] + 1
  #tmp = vfCrossCorr[zz,:]
  #maximaLoc = np.r_[True, tmp[1:] > tmp[:-1]] & np.r_[tmp[:-1] > tmp[1:], True]
  #maximaLoc = np.squeeze(np.argwhere(maximaLoc).T)
  # maximum values
  maxima = vfCrossCorr[zz,maximaLoc]

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

tauDistance = 0.01/tstar  # to make sure maximas are near each other
for zz in np.arange(0,nz-1):
  tauFirstNext = vfFirstMaxima[zz+1,0]
  tauSecondNext = vfSecondMaxima[zz+1,0]
  tauThirdNext = vfThirdMaxima[zz+1, 0]

  # Use first max if it is close
  if np.abs(tau[zz] - tauFirstNext) < tauDistance:
    tau[zz+1] = tauFirstNext
    zeta[zz+1] = vfFirstMaxima[zz+1,1]
  # if it's too far, try the second max
  elif np.abs(tau[zz] - tauSecondNext) < tauDistance:
    tau[zz+1] = tauSecondNext
    zeta[zz+1] = vfSecondMaxima[zz+1,1]
  # if it's too far, try the third
  elif np.abs(tau[zz] - tauThirdNext) < tauDistance:
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
vfFig = plt.figure(figsize=(3.25,1.625))
plt.imshow(vfCrossCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]],
  vmin=-1., vmax=1., cmap='coolwarm')
#plt.colorbar()

#plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'w-', alpha=0.4)
#plt.plot(vfSecondMaxima[:,0], vfSecondMaxima[:,1], 'r-', alpha=0.4)
#plt.plot(vfThirdMaxima[:,0], vfThirdMaxima[:,1], 'g-', alpha=0.4)
plt.plot(tau, zeta, 'k-')

plt.plot(tau, yFit, 'w--')
#cTxtString = r"$c_* = \Delta_z = %.3f\Delta t$" % p
cTxtString = r"$c_* = %.3f$" % p
plt.text(2, 10, cTxtString, fontsize=12)

#plt.xlabel(r'$\Delta t_*$')
plt.xlabel(r'$\nu \Delta t / (2a)^2$')
plt.xlim([0, 6])
#plt.ylabel(r'$\Delta z_*$', rotation=0)
plt.ylabel(r'$\Delta z/2a$', rotation=90)
plt.ylim([dz[0], dz[-1]])
plt.gca().yaxis.set_label_coords(-0.15, 0.45)

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(5))
plt.gca().yaxis.set_minor_locator(MultipleLocator(2.5))

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
