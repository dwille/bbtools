#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                Surface Plotting"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

## VOLUME FRACTION ##
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T

# Plot
vFracFig = plt.figure()

minVal = np.floor(100*np.amin(vFrac))/100
maxVal = np.ceil(100*np.amax(vFrac))/100

plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=minVal, vmax=maxVal)

cbar = plt.colorbar()
plt.title(r'$Volume\ Fraction$')
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "volume-fraction"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# ## NUMBER DENSITY ##
#nDensFile = datadir + "number-density"
#numDens = np.genfromtxt(nDensFile).T
# nDensFig = plt.figure()
# plt.imshow(numDens, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], evalZ[0], evalZ[-1]])
# plt.colorbar(format="%.1e")
# plt.title(r"$n$")
# plt.xlabel(r"$t - t_0$")
# plt.xticks(np.floor(np.arange(time[0], time[-1], 1000)))
# plt.ylabel(r'$z/a$')
# 
# imgname = imgdir + "number-density"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
# #plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## u_PART ##
upFig = plt.figure()

upFile = datadir + "part-u"
up = np.genfromtxt(upFile).T

minVal = np.floor(100*np.amin(up))/100
maxVal = np.ceil(100*np.amax(up))/100
maxVal = np.max((np.abs(minVal), np.abs(maxVal)))

plt.imshow(up, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-maxVal, vmax=maxVal)

cbar = plt.colorbar()
plt.title(r"$U_p\ [mm/ms]$")
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "part-u"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## v_PART ##
vpFig = plt.figure()

vpFile = datadir + "part-v"
vp = np.genfromtxt(vpFile).T

minVal = np.floor(100*np.amin(vp))/100
maxVal = np.ceil(100*np.amax(vp))/100
maxVal = np.max((np.abs(minVal), np.abs(maxVal)))

plt.imshow(vp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-maxVal, vmax=maxVal)

cbar = plt.colorbar()
plt.title(r"$V_p\ [mm/ms]$")
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "part-v"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## W_PART ##
wpFig = plt.figure()

wpFile = datadir + "part-w"
wp = np.genfromtxt(wpFile).T

minVal = np.floor(100*np.amin(wp))/100
maxVal = np.ceil(100*np.amax(wp))/100
maxVal = np.max((np.abs(minVal), np.abs(maxVal)))

plt.imshow(wp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-maxVal, vmax=maxVal)

cbar = plt.colorbar()
plt.title(r"$W_p\ [mm/ms]$")
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "part-w"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Kinetic Energy ##
keFig = plt.figure()

keFile = datadir + "part-kinetic-energy"
ke = np.genfromtxt(keFile).T

maxVal = np.ceil(1000*np.amax(ke))/1000
maxVal = np.max((np.abs(minVal), np.abs(maxVal)))
maxVal = np.amax(ke)

plt.imshow(ke, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=0, vmax=maxVal)

cbar = plt.colorbar()
plt.title(r"$\kappa\ [mm^2/ms^2]$")
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "part-ke"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
