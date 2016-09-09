#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Avg Waveform"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = np.mean(np.diff(evalZ))
dt = np.mean(np.diff(time))
halfNz = int(np.floor(0.5*nz))

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# loop params
zs = 0
#zloop = halfNz
zloop = nz

# Set up storage for xcorr: (dtau_index, dz_index)
vfFirstMaxima  = np.zeros((zloop,2))
vfSecondMaxima = np.zeros((zloop,2))
vfThirdMaxima  = np.zeros((zloop,2))
vfFourthMaxima  = np.zeros((zloop,2))

# Cross correlate slices at constant z -- over zloop 
for zz in np.arange(0,zloop):      # distance away from start
  # correctly loop through domain
  if zs + zz >= nz:
    zInd = zs + zz - nz
  else:
    zInd = zs + zz

  # length of result is ceil(length(time)/2)
  vfCrossCorr = CrossCorrelationFFT(vFrac[zs,:], vFrac[zInd,:])
  if (zz == 0):
    norm = vfCrossCorr[0]
  vfCrossCorr /= norm

  # Find maxima locations by using where second derivative changes
  maximaLoc = (np.diff(np.sign(np.diff(vfCrossCorr))) < 0).nonzero()[0] + 1
  # maximum values
  maxima = vfCrossCorr[maximaLoc]

  # sort the maxima
  if np.size(maximaLoc) == 0:       # if no maxima, give it trash
    vfFirstMaxima[zz,0] = np.nan
    vfFirstMaxima[zz,1] = zz
    vfSecondMaxima[zz,0] = np.nan
    vfSecondMaxima[zz,1] = zz
    vfThirdMaxima[zz,0] = np.nan
    vfThirdMaxima[zz,1] = zz
    vfFourthMaxima[zz,0] = np.nan
    vfFourthMaxima[zz,1] = zz
  elif np.size(maximaLoc) > 0:      # If at least one, set the first
    vfFirstMaxima[zz,0] = maximaLoc[0]
    vfFirstMaxima[zz,1] = zz

    if np.size(maximaLoc) > 1:      # If at least two, set the second
      vfSecondMaxima[zz,0] = maximaLoc[1]
      vfSecondMaxima[zz,1] = zz

      if np.size(maximaLoc) > 2:    # If at least 3, set the third
        vfThirdMaxima[zz,0] = maximaLoc[2]
        vfThirdMaxima[zz,1] = zz

        if np.size(maximaLoc) > 3:  # If at least 4, set 4th
          vfFourthMaxima[zz,0] = maximaLoc[3]
          vfFourthMaxima[zz,1] = zz

  if vfCrossCorr[0] > vfCrossCorr[1]:  # If [0] is a max, shift them all
    if np.size(maximaLoc) > 2:
      vfFourthMaxima[zz,0] = vfThirdMaxima[zz,0]
      vfFourthMaxima[zz,1] = vfThirdMaxima[zz,1]

    if np.size(maximaLoc) > 2:
      vfThirdMaxima[zz,0] = vfSecondMaxima[zz,0]
      vfThirdMaxima[zz,1] = vfSecondMaxima[zz,1]

    if np.size(maximaLoc) > 1:
      vfSecondMaxima[zz,0] = vfFirstMaxima[zz,0]
      vfSecondMaxima[zz,1] = vfFirstMaxima[zz,1]

    vfFirstMaxima[zz,0] = 0
    vfFirstMaxima[zz,1] = zz

# Find actual maxima that traces the waveform      
tauInd = -np.ones(zloop)
zetaInd = -np.ones(zloop)
tauInd[0] = vfFirstMaxima[0,0]
zetaInd[0] = vfFirstMaxima[0,1]
timeIndComp = 0.01/dt   # "close enough"
for zz in np.arange(0, zloop-1):
  tauIndFirstNext  =  vfFirstMaxima[zz+1,0]
  tauIndSecondNext = vfSecondMaxima[zz+1,0]
  tauIndThirdNext  =  vfThirdMaxima[zz+1,0]
  tauIndFourthNext  =  vfFourthMaxima[zz+1,0]

  # Use first max if it is close enough
  if np.abs(tauInd[zz] - tauIndFirstNext) < timeIndComp:
    tauInd[zz+1] = tauIndFirstNext
    zetaInd[zz+1] = vfFirstMaxima[zz+1,1]
  # if it's too far, try the second max
  elif np.abs(tauInd[zz] - tauIndSecondNext) < timeIndComp:
    tauInd[zz+1] = tauIndSecondNext
    zetaInd[zz+1] = vfSecondMaxima[zz+1,1]
  # if it's too far, try the third
  elif np.abs(tauInd[zz] - tauIndThirdNext) < timeIndComp:
    tauInd[zz+1] = tauIndThirdNext
    zetaInd[zz+1] = vfThirdMaxima[zz+1,1]
  # .... and the fourth
  elif np.abs(tauInd[zz] - tauIndFourthNext) < timeIndComp:
    tauInd[zz+1] = tauIndFourthNext
    zetaInd[zz+1] = vfFourthMaxima[zz+1,1]
  # if that is not good, quit
  else:
    break;

tauInd = tauInd[tauInd > -1]
nmax = np.max(tauInd)
zetaInd = zetaInd[zetaInd > -1]

# Average waveforms together
ntmax = nt - nmax
avgFormTime = np.zeros(ntmax)
for zz in np.arange(0,np.size(tauInd)):
  # correctly loop through domain
  if zs + zz >= nz:
    zInd = zs + zz - nz
  else:
    zInd = zs + zz

  avgFormTime += vFrac[zInd,tauInd[zz]:(tauInd[zz]+ntmax)]/np.size(tauInd)

avgFormSpace = np.zeros(nz)
for zz in np.arange(0,np.size(tauInd)):
  # correctly loop through domain
  if zs + zz >= nz:
    zInd = zs + zz - nz
  else:
    zInd = zs + zz

  ind = np.arange(nz)
  ind = np.roll(ind, int(-zetaInd[zz]))
  ind = np.roll(ind, halfNz)

  avgFormSpace += vFrac[ind,tauInd[zz]]/np.size(tauInd)

# Plot
fig1 = plt.figure(figsize=(4,6))
subs = 310

ax3 = fig1.add_subplot(subs+3)
plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[0, np.size(time), 0, nz])
plt.plot(vfFirstMaxima[:,0], vfFirstMaxima[:,1], '.')
plt.plot(vfSecondMaxima[:,0], vfSecondMaxima[:,1], '.')
plt.plot(vfThirdMaxima[:,0], vfThirdMaxima[:,1], '.')
plt.plot(vfFourthMaxima[:,0], vfFourthMaxima[:,1], '.')
plt.plot(tauInd[tauInd > -1], zetaInd[zetaInd > -1], 'k.')
ax3.set_xlim([0, 500])

ax1 = fig1.add_subplot(subs+1)
plt.plot(time[0:ntmax],avgFormTime)
plt.xlim([0,4])

ax2 = fig1.add_subplot(subs+2)
plt.plot(evalZ, avgFormSpace)#, evalZ)
ax2.set_xlim([evalZ[0], evalZ[-1]])

imgname = imgdir + "avg-waveform"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
