#!/usr/bin/env python2
from setup import *
os.system('clear')

from skimage.io import imread
from skimage import data_dir
from skimage.transform import radon, rescale

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
vFrac = vFrac.T

sinogram = radon(vFrac)

## Save Data to File
savefile = datadir + "radon-xform"
with open(savefile, 'wb') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerows(sinogram)

print "\n      ...Done!"
