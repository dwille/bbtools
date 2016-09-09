#!/usr/bin/env python2
from setup import *
os.system('clear')
from matplotlib.ticker import MultipleLocator

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                Surface Plotting"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])
nu = .01715   ## mm^2/ms XXX

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

## VOLUME FRACTION ##
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# get avg volume fraction
domX = 42 # XXX
domY = 42 # XXX
domZ = 126 # XXX
avgVolumeFraction = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)

# get phase averaged data
phaseVelDir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVelDir, skip_header=1)

for nn in np.arange(0,np.size(phaseVelData, 0)):
  currN = phaseVelData[nn,0]
  currRho = phaseVelData[nn,1]
  if (nparts == currN) & (rho == currRho):
    phaseVel = phaseVelData[nn,2]
phaseVel *= 1000  # mm / ms -> mm / s

# Plot
vFracFig = plt.figure(figsize=(3.25,1.625))

minVal = np.floor(100*np.amin(vFrac))/100
maxVal = np.ceil(100*np.amax(vFrac))/100
upperDiff = maxVal - avgVolumeFraction
lowerDiff = avgVolumeFraction - minVal
maxDiff = np.max([upperDiff, lowerDiff])
maxDiff = np.floor(maxDiff * 100.) / 100.

# Highlight Fluid
#vFrac[vFrac >= avgVolumeFraction] = avgVolumeFraction

# Highight Solid
#vFrac[vFrac <= avgVolumeFraction] = avgVolumeFraction

#tau = (2.*partR) / phaseVel       ## mm / (mm/s) = s
# Non-dimensionalize
evalZ /= (2.*partR)               ## mm / mm
tau = (2.*partR)*(2.*partR)/nu    ## mm^2/(mm^2/ms) = ms
tau /= 1000                       ## ms -> s
time /= tau                       ## s / s

plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=avgVolumeFraction - maxDiff, vmax=avgVolumeFraction + maxDiff,
  cmap='coolwarm')

cbar = plt.colorbar()
plt.xlabel(r"$\nu t/(2a)^2$")
plt.ylabel(r'$z/2a$', rotation=90)
plt.gca().yaxis.set_label_coords(-0.15, 0.5)

xEnd = time[-1]
xEnd = 10.*np.floor(xEnd / 10.)
xEnd = 14 # XXX
plt.xlim([0, xEnd])
plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(5))
plt.gca().yaxis.set_minor_locator(MultipleLocator(2.5))

imgname = imgdir + "volume-fraction"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')

print "I'm lazy, let's exit here. (line 66)"
sys.exit()

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
 
## Kinetic Energy ##
keFig = plt.figure(figsize=(3.25,1.625))

keFile = datadir + "part-kinetic-energy"
ke = np.genfromtxt(keFile).T

#ke = np.multiply(ke, vFrac)

maxVal = np.ceil(1000*np.amax(ke))/1000
maxVal = np.max((np.abs(minVal), np.abs(maxVal)))
maxVal = np.amax(ke)

plt.imshow(ke, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=0, vmax=maxVal, cmap="viridis")

cbar = plt.colorbar()
#plt.title(r"$\kappa\ [mm^2/ms^2]$")
plt.xlabel(r"$t\ [s]$")
plt.ylabel(r'$z\ [mm]$')

xEnd = time[-1]
plt.xlim([0, xEnd])
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))

imgname = imgdir + "part-ke"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
