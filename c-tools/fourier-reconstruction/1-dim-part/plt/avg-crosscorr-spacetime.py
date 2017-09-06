#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Avg Crosscorrelate"
print ""

# Setup simulation parameters
(partR,_,_,_, simdir, tstart) = simParams(sys)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

## Crosscorrelation of volume fraction ##
# crosscorrelate phi(z_i), phi(z_j), e.g.
# phi(z_i), phi(z_i + dz)

# each of (0,nz,nt) is a correlation at a starting zs
# we will average these over (nz,:,:)
# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]  # ignore time before tstart

vfCrossCorr = np.zeros((nz,nt))
temp = np.zeros((nz,nt))
for zs in np.arange(0,nz):          # starting index, z_i
  for zz, zval in enumerate(evalZ): # correlation index, z_j
    # correctly loop through domain
    if zs + zz >= nz:
      zInd = zs + zz - nz
    else:
      zInd = zs + zz

    # length of result is ceil(length(time)/2)
    temp[zz,:] = CrossCorrelationFFT(vFrac[zs,:], vFrac[zInd,:])

  temp /= temp[0,0]
  vfCrossCorr += temp/nz

## Save data to file
savefile = datadir + "z-averaged-vfrac-xcorr"
with open(savefile, 'wb') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerows(vfCrossCorr)

# ## Crosscorrelation of wp ##
# # each of (0,nz,nt) is a correlation at a starting zs
# # we will average these over (nz,:,:)
# wp = np.genfromtxt(datadir + "part-w").T[:,tsInd:]  # ignore time before tstart
# 
# wpCrossCorr = np.zeros((nz,nt))
# temp = np.zeros((nz,nt))
# for zs in np.arange(0,nz):
#   for zz, zval in enumerate(evalZ):
#     # correctly loop through domain
#     if zs + zz >= nz:
#       zInd = zs + zz - nz
#     else:
#       zInd = zs + zz
# 
#     # length of result is ceil(length(time)/2)
#     temp[zz,:] = CrossCorrelationFFT(wp[zs,:], wp[zInd,:])
# 
#   temp /= temp[0,0]
#   wpCrossCorr += temp/nz
# 
# 
# ## Save data to file
# savefile = datadir + "z-averaged-wp-xcorr"
# with open(savefile, 'wb') as outfile:
#   out = csv.writer(outfile, delimiter=' ')
#   out.writerows(wpCrossCorr)
print "      Done!"
