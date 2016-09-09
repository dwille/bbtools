#!/usr/bin/env python2
from setup import *
os.system('clear')

print 
print " ---- 1D Fourier Reconstruction ---- "
print "          Fluct Magnitude"
print

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)
nparts = int(simdir.partition('/')[0])
domX = 42
domY = 42
domZ = 126
avgVolumeFraction = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Find max / mins
vfMax_1 = np.zeros((nz,3))
vfMax_2 = np.zeros((nz,3))
vfMax_3 = np.zeros((nz,3))
vfMin_1 = np.zeros((nz,3))
vfMin_2 = np.zeros((nz,3))
vfMin_3 = np.zeros((nz,3))
for zz in np.arange(0,nz):
  # Find max/min locs by finding where second derivative changes
  # Does NOT pick up index 0, if it is a max
  secondDerivative = np.diff(np.sign(np.diff(vFrac[zz,:])))
  maxLoc = (secondDerivative < 0).nonzero()[0] + 1
  minLoc = (secondDerivative > 0).nonzero()[0] + 1

  # Max/min values
  maxima = vFrac[zz,maxLoc]
  minima = vFrac[zz,minLoc]

  # Max/min indices
  maxInd = np.argmax(maxima)
  minInd = np.argmin(minima)

  # Sort the maxima
  if np.size(maxLoc) == 0:       # If no maxima, give it trash data
    vfMax_1[zz,0] = np.nan
    vfMax_1[zz,1] = evalZ[zz]
    vfMax_1[zz,2] = np.nan

    vfMax_2[zz,0] = np.nan
    vfMax_2[zz,1] = evalZ[zz]
    vfMax_2[zz,2] = np.nan
  
  elif np.size(maxLoc) > 0:      # If at least one, set the first
    vfMax_1[zz,0] = time[maxLoc[0]]
    vfMax_1[zz,1] = evalZ[zz]
    vfMax_1[zz,2] = vFrac[zz,maxLoc[0]]

    if np.size(maxLoc) > 1:      # If at least two, set the second
      vfMax_2[zz,0] = time[maxLoc[1]]
      vfMax_2[zz,1] = evalZ[zz]
      vfMax_2[zz,2] = vFrac[zz,maxLoc[1]]

      if np.size(maxLoc) > 2:    # If at least 3, set the third
        vfMax_3[zz,0] = time[maxLoc[2]]
        vfMax_3[zz,1] = evalZ[zz]
        vfMax_3[zz,2] = vFrac[zz,maxLoc[2]]

  if vFrac[zz,0] > vFrac[zz,1]:   # If [0] is a max, shift them all
    if np.size(maxLoc) > 2:
      vfMax_3[zz,0] = vfMax_2[zz,0]
      vfMax_3[zz,1] = vfMax_2[zz,1]
      vfMax_3[zz,2] = vfMax_2[zz,2]

    if np.size(maxLoc) > 1:
      vfMax_2[zz,0] = vfMax_1[zz,0]
      vfMax_2[zz,1] = vfMax_1[zz,1]
      vfMax_2[zz,2] = vfMax_1[zz,2]

    vfMax_1[zz,0] = time[0]
    vfMax_1[zz,1] = evalZ[zz]
    vfMax_1[zz,2] = vFrac[zz,0]


fig1 = plt.figure()
# Plotting

minVal = np.floor(100*np.amin(vFrac))/100
maxVal = np.ceil(100*np.amax(vFrac))/100
upperDiff = maxVal - avgVolumeFraction
lowerDiff = avgVolumeFraction - minVal
maxDiff = np.max([upperDiff, lowerDiff])

# vfrac color
ax1 = fig1.add_subplot(211)
plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=avgVolumeFraction + maxDiff, vmax=avgVolumeFraction - maxDiff,
  cmap='seismic')
plt.colorbar()
ax1.set_xlim([0,2])
ax1.set_ylim([-63,63])

# Maxima
ax1.plot(vfMax_1[:,0], vfMax_1[:,1], 'k.')
ax1.plot(vfMax_2[:,0], vfMax_2[:,1], 'r.')
ax1.plot(vfMax_3[:,0], vfMax_3[:,1], 'g.')

# Hist
ax2 = fig1.add_subplot(212)

nBins = 25
hist, bin_edges = np.histogram(np.reshape(vFrac,(np.size(vFrac),-1)), nBins)
dBin = np.mean(np.diff(bin_edges))
hist = hist/(np.sum(hist) * dBin)      # convert from counts to pdf
bin_centers = bin_edges[0:-1] + dBin
ax2.plot(bin_centers, hist)
 
# Find confidence intervals
cdf = np.cumsum(hist)*dBin
ind_01 = (cdf < 0.01).nonzero()
ind_05 = (cdf < 0.05).nonzero()
ind_95 = (cdf > 0.95).nonzero()
ind_99 = (cdf > 0.99).nonzero()

# stats
print "Min  = %.3f" % np.min(vFrac)
print ".01  = %.3f" % bin_centers[ind_01][-1]
print ".05  = %.3f" % bin_centers[ind_05][-1]
print "-Std = %.3f" % (avgVolumeFraction - np.std(vFrac))
print "Mean = %.3f" % avgVolumeFraction
print "sdev = %.3f" % np.std(vFrac)
print "+Std = %.3f" % (avgVolumeFraction + np.std(vFrac))
print ".95  = %.3f" % bin_centers[ind_95][0]
print ".99  = %.3f" % bin_centers[ind_99][0]
print "Max  = %.3f" % np.max(vFrac)


# Save
imgname = imgdir + "vfrac-flucts"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
print "      Done!"
