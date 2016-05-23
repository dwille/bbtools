#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "             Autocorrelation in Time"
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

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Autocorrelation of volume fraction
vfFig = plt.figure()
vfAutoCorr = np.zeros((nz, nt))
vfMaxima = np.zeros((nt,3))
vfFirstMaxima = np.zeros((nt,3))
for tt,tval in enumerate(time):

  # length of result is ceil(length(time)/2)
  vfAutoCorr[:,tt] = AutoCorrelationFFT(vFrac[:,tt])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(vfAutoCorr[:,tt]))) < 0).nonzero()[0] + 1
  maxima = vfAutoCorr[maximaLoc,tt]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfMaxima[tt,0] = time[tt]
    vfMaxima[tt,1] = np.nan
    vfMaxima[tt,2] = np.nan
  else:
    vfMaxima[tt,0] = time[tt]
    vfMaxima[tt,1] = dz[maximaLoc[maxInd]]
    vfMaxima[tt,2] = vfAutoCorr[maximaLoc[maxInd],tt]

    vfFirstMaxima[tt,0] = time[tt]
    vfFirstMaxima[tt,1] = dz[maximaLoc[0]]
    vfFirstMaxima[tt,2] = vfAutoCorr[maximaLoc[0],tt]

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])
plt.colorbar()
plt.plot(vfMaxima[:,0], vfMaxima[:,1], '--', color='0.75')
plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'k--')

plt.xlabel(r"$t - t_0\ [s]$")
plt.ylabel(r"$\Delta z\ [mm]$",rotation=0)
plt.title(r"$\langle \alpha(z) \alpha(z + \Delta z) \rangle $")
plt.xlim([time[0], time[-1]])
plt.ylim([dz[0], dz[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

imgname = imgdir + "autocorr-space-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# # Autocorrelation of particle velocity - space
# wpFig = plt.figure()
# wpAutoCorr = np.zeros((nz,nt))
# wpMaxima = np.zeros((nz,3))
# wpFirstMaxima = np.zeros((nt,3))
# for tt,tval in enumerate(time):
# 
#   # length of result is ceil(length(time)/2)
#   wpAutoCorr[:,tt] = AutoCorrelationFFT(wp[:,tt])
# 
#   # find maxima
# #  maximaLoc = (np.diff(np.sign(np.diff(wpAutoCorr[:,tt]))) < 0).nonzero()[0] + 1
# #  maxima = wpAutoCorr[maximaLoc,tt]
# #  maxInd = np.argmax(maxima)
# #  if np.size(maximaLoc) == 0:
# #    wpMaxima[tt,0] = time[tt]
# #    wpMaxima[tt,1] = np.nan
# #    wpMaxima[tt,2] = np.nan
# #  else:
# #    wpMaxima[tt,0] = time[tt]
# #    wpMaxima[tt,1] = dz[maximaLoc[maxInd]]
# #    wpMaxima[tt,2] = wpAutoCorr[maximaLoc[maxInd],tt]
# #
# #    wpFirstMaxima[tt,0] = time[tt]
# #    wpFirstMaxima[tt,1] = dz[maximaLoc[0]]
# #    wpFirstMaxima[tt,2] = wpAutoCorr[maximaLoc[0],tt]
# 
# plt.imshow(wpAutoCorr, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], dz[0], dz[-1]])
# plt.colorbar()
# 
# #plt.plot(wpMaxima[:,0], wpMaxima[:,1], 'ko--')
# 
# plt.xlabel(r"$t - t_0\ [s]$")
# plt.ylabel(r"$\Delta z\ [mm]$",rotation=0)
# plt.title(r"$\langle w_p(z) w_p(z + \Delta z) \rangle $")
# plt.xlim([time[0], time[-1]])
# plt.ylim([dz[0], dz[-1]])
# plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
# 
# imgname = imgdir + "autocorr-space-wp"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
print "      Done!"
