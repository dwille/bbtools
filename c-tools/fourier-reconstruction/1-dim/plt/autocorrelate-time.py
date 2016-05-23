#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Autocorrelation in Space"
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

# Autocorrelation of volume fraction
vfFig = plt.figure()
vfAutoCorr = np.zeros((nz, nt))
vfMaxima = np.zeros((nz,3))
vfFirstMaxima = np.zeros((nz,3))
for zz, zval in enumerate(evalZ):

  # length of result is ceil(length(time)/2)
  vfAutoCorr[zz,:] = AutoCorrelationFFT(vFrac[zz,:])

  # find maxima
  maximaLoc = (np.diff(np.sign(np.diff(vfAutoCorr[zz,:]))) < 0).nonzero()[0] + 1
  maxima = vfAutoCorr[zz,maximaLoc]
  maxInd = np.argmax(maxima)
  if np.size(maximaLoc) == 0:
    vfMaxima[zz,0] = np.nan
    vfMaxima[zz,1] = evalZ[zz]
    vfMaxima[zz,2] = np.nan
  else:
    vfMaxima[zz,0] = time[maximaLoc[maxInd]]
    vfMaxima[zz,1] = evalZ[zz]
    vfMaxima[zz,2] = vfAutoCorr[zz,maximaLoc[maxInd]]

    vfFirstMaxima[zz,0] = time[maximaLoc[0]]
    vfFirstMaxima[zz,1] = evalZ[zz]
    vfFirstMaxima[zz,2] = vfAutoCorr[zz,maximaLoc[0]]

tauVal = np.mean(vfFirstMaxima[:,0])
tauMean = np.mean(vfFirstMaxima[:,2])
print ""
print "      Mean = %.2f at tau = %.4f" % (tauMean, tauVal)
print "      omega = 2pi/tau = %.4f" % (2.*np.pi/tauVal)
print "      Freq = 1/tau = %.4f" % (1./tauVal)

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
plt.colorbar()
#plt.plot(vfMaxima[:,0], vfMaxima[:,1], '--', color='0.75')
plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], 'k--')

plt.plot([tauVal, tauVal], [evalZ[0], evalZ[-1]], 'k:')
txtString = r"$\bar \tau = %.4f$" % tauVal
plt.text(1, 0, txtString, fontsize=12)

plt.xlabel(r"$\tau\ [s]$")
plt.ylabel(r"$z\ [mm]$",rotation=0)
plt.title(r"$\langle \alpha(t) \alpha(t + \tau) \rangle $")
plt.xlim([time[0], time[-1]])
plt.ylim([evalZ[0], evalZ[-1]])
plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))

imgname = imgdir + "autocorr-time-vf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# # Autocorrelation of wp -- time
# wpFile = datadir + "part-w"
# wp = np.genfromtxt(wpFile).T[:,tsInd:]
#
# wpFig = plt.figure()
# wpAutoCorr = np.zeros((nz, nt))
# wpMaxima = np.zeros((nz,3))
# for zz, zval in enumerate(evalZ):
#   # length of result is ceil(length(time)/2)
#   wpAutoCorr[zz,:] = AutoCorrelationFFT(wp[zz,:])
# 
#   # find maxima
#   maximaLoc = (np.diff(np.sign(np.diff(wpAutoCorr[zz,:]))) < 0).nonzero()[0] + 1
#   if np.size(maximaLoc) == 0:
#     wpMaxima[zz,0] = np.nan
#     wpMaxima[zz,1] = evalZ[zz]
#     wpMaxima[zz,2] = np.nan
#   else:
#     wpMaxima[zz,0] = time[maximaLoc[0]]
#     wpMaxima[zz,1] = evalZ[zz]
#     wpMaxima[zz,2] = wpAutoCorr[zz,maximaLoc[0]]
# 
# plt.imshow(wpAutoCorr, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
# plt.colorbar()
# plt.plot(wpMaxima[:,0], wpMaxima[:,1], 'k--')
# 
# plt.xlabel(r"$\tau\ [s]$")
# plt.ylabel(r"$z\ [mm]$",rotation=0)
# plt.title(r"$\langle w_p(t) w_p(t + \tau) \rangle $")
# plt.xlim([time[0], time[-1]])
# plt.ylim([evalZ[0], evalZ[-1]])
# plt.xticks(np.floor(np.arange(time[0], time[-1], 1)))
# 
# imgname = imgdir + "autocorr-time-wp"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
# 
print "\n      ...Done!"
