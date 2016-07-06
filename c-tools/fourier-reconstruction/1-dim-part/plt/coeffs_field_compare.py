#!/usr/bin/env python2
from setup import *
os.system('clear')

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                  Field Comparison"
print ""

# SIMULATION PARAMETERS
domX = 42
domY = 42
domZ = 126
xn = 160
yn = 160
zn = 480
z = np.linspace(-0.5*domZ, 0.5*domZ, zn)
z = z[np.newaxis]

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

nparts = float(simdir.partition('/')[0])

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# set up output file paths
infoFile = datadir + "info"
n15File = datadir + "volume-fraction-15"
n15 = np.genfromtxt(n15File).T
# norm of max coeffs
norm15 = np.sqrt(np.sum(n15**2., 0))

# pull the rest
coeff_list = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
  "11", "12", "13", "14"]
err = np.zeros(np.size(coeff_list) + 1)
for pp,pull in enumerate(coeff_list):
  currFile = datadir + "volume-fraction-" + pull
  try: 
    nCurr = np.genfromtxt(currFile).T
    # take difference with full reconstruction
    diff = np.sqrt(np.sum((n15 - nCurr)**2., 0))
    # normalize and take mean over time
    err[pp] = np.mean(diff / norm15)
  except:
    print "File " + currFile + " not found."


fig1 = plt.figure()
plt.plot(np.arange(0, 16), err, 'ko')

plt.xlabel("Orders Used")
plt.ylabel("Time-Averaged RMS with 15")
plt.xlim([-0.5, 15.5])
plt.ylim([0., 0.2])


imgname = imgdir + "vf-field-compare"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
