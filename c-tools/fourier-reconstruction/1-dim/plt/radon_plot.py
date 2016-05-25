#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "             Radon Xform Plot"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
#(root, simdir, datadir, imgdir) = directoryStructureDevel(simdir)
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = np.mean(np.diff(evalZ))
dt = np.mean(np.diff(time))
dzdt = dz/dt

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Pull output data -- each column is a different time
savefile = datadir + "radon-xform"
radonXForm = np.genfromtxt(savefile)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.imshow(radonXForm, extent=(0, 180, 0, radonXForm.shape[0]), aspect='auto')

ax1.set_xticks([0,36,72,108,144,180])
currTicks = ax1.xaxis.get_ticklocs()
scaledTicks = np.tan(np.deg2rad(currTicks)) * dzdt
scaledTicks = np.round(scaledTicks*100.)/100
ax1.set_xticklabels(scaledTicks)

ax1.set_xlabel('Projection angle (deg)')
ax1.set_ylabel('Projection position (pixels)')

# this is obnoxious. need to scale plot to be non-monotonic / wraparound
c = 119.09
cloc = np.rad2deg(np.arctan(119.09/dzdt))
ax1.plot([cloc, cloc], [0, radonXForm.shape[0]], 'k--')



#plt.xlim([45,90])

imgname = imgdir + "radon-xform"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
