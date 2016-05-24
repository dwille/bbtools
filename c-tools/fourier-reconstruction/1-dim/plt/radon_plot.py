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
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Pull output data -- each column is a different time
savefile = datadir + "radon-xform"
radonXForm = np.genfromtxt(savefile)

fig1 = plt.figure()
plt.imshow(radonXForm, extent=(0, 180, 0, radonXForm.shape[0]), aspect='auto')
plt.xlabel('Projection angle (deg)')
plt.ylabel('Projection position (pixels)')

imgname = imgdir + "radon-xform"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
