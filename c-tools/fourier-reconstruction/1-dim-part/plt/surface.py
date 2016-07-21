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
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# get avg volume fraction
domX = 42
domY = 42
domZ = 126
nparts = int(simdir.partition('/')[0])
avgVolumeFraction = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)

# Plot
vFracFig = plt.figure(figsize=(6,3))

minVal = np.floor(100*np.amin(vFrac))/100
maxVal = np.ceil(100*np.amax(vFrac))/100
upperDiff = maxVal - avgVolumeFraction
lowerDiff = avgVolumeFraction - minVal
maxDiff = np.max([upperDiff, lowerDiff])

# Highlight Fluid
vFrac[vFrac >= avgVolumeFraction] = avgVolumeFraction

# Highight Solid
#vFrac[vFrac <= avgVolumeFraction] = avgVolumeFraction

plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=avgVolumeFraction + maxDiff, vmax=avgVolumeFraction - maxDiff,
  cmap='seismic')

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
#numDens = np.genfromtxt(nDensFile).T[:,tsInd:]
# nDensFig = plt.figure(figsize=(6,3))
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

# ## u_PART ##
upFig = plt.figure(figsize=(6,3))

upFile = datadir + "part-u"
up = np.genfromtxt(upFile).T[:,tsInd:]
vpFile = datadir + "part-v"
vp = np.genfromtxt(vpFile).T[:,tsInd:]

minValU = np.floor(100*np.amin(up))/100
maxValU = np.ceil(100*np.amax(up))/100
minValV = np.floor(100*np.amin(vp))/100
maxValV = np.ceil(100*np.amax(vp))/100
maxVal = np.max((np.abs(minValU), np.abs(maxValU), np.abs(maxValV), np.abs(minValV)))

plt.imshow(up, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-maxVal, vmax=maxVal, cmap="seismic")

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
vpFig = plt.figure(figsize=(6,3))

plt.imshow(vp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-maxVal, vmax=maxVal, cmap="seismic")

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

## u perp ##
uperp = np.sqrt(up**2 + vp**2)
uperp -= np.mean(uperp)

uperpFig = plt.figure(figsize=(6,3))

plt.imshow(uperp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]], 
  vmin=-np.max(np.abs(uperp)), vmax=np.max(np.abs(uperp)), cmap="plasma")

cbar = plt.colorbar()
plt.title(r"$u_\perp^\prime\ [mm/ms]$")
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "part-uperp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## W_PART ##
wpFig = plt.figure(figsize=(6,3))

wpFile = datadir + "part-w"
wp = np.genfromtxt(wpFile).T[:,tsInd:]

minVal = np.floor(100*np.amin(wp))/100
maxVal = np.ceil(100*np.amax(wp))/100
maxVal = np.max((np.abs(minVal), np.abs(maxVal)))

plt.imshow(wp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-maxVal, vmax=maxVal, cmap='seismic')

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
# 
# ## Kinetic Energy ##
# keFig = plt.figure()
# 
# keFile = datadir + "part-kinetic-energy"
# ke = np.genfromtxt(keFile).T
# 
# maxVal = np.ceil(1000*np.amax(ke))/1000
# maxVal = np.max((np.abs(minVal), np.abs(maxVal)))
# maxVal = np.amax(ke)
# 
# plt.imshow(ke, origin="lower", aspect="auto", interpolation="none",
#   extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
#   vmin=0, vmax=maxVal)
# 
# cbar = plt.colorbar()
# plt.title(r"$\kappa\ [mm^2/ms^2]$")
# plt.xlabel(r"$t\ [s]$")
# plt.ylabel(r'$z\ [mm]$')
# 
# xEnd = time[-1]
# plt.xlim([0, xEnd])
# plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))
# 
# imgname = imgdir + "part-ke"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
# #plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
