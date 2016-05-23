#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                Slices Plotting"
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

##
 # Read Volume Fraction
 ##
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

tt = 250
zz = np.floor(np.size(evalZ)/2)
labelx = -0.15

##
 # Plot one volume fraction slice
 ##
vFracFig = plt.figure(figsize=(6,4))
gs = matplotlib.gridspec.GridSpec(2,3)

# fixed time slice
ax1 = vFracFig.add_subplot(gs[0,0])
ax1.plot(vFrac[:,tt], evalZ, 'b', linewidth=2)

ax1.set_ylim([np.min(evalZ), np.max(evalZ)])
ax1.set_xticks([0.3, 0.35, 0.40])
ax1.set_xlabel(r'$\phi$', fontsize=14)
ax1.set_ylabel(r'$z\ [mm]$', fontsize=14)

# fixed location slice
ax2 = vFracFig.add_subplot(gs[1,1:3])
plt.plot(time, vFrac[zz,:], 'k', linewidth=2)

xEnd = time[-1]
plt.xlim([0, xEnd])
ax2.set_xlabel(r"$t\ [s]$", fontsize=14)
ax2.set_ylabel(r"$\phi$", fontsize=14)
plt.xticks(np.floor(np.arange(0, xEnd+0.01, 1)))
ax2.set_yticks([0.3, 0.35, 0.40])

# contour plot
ax3 = vFracFig.add_subplot(gs[0,1:3])

im = ax3.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
      extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
      vmin=0.3, vmax=0.4)

ax3.plot([time[tt], time[tt]], [evalZ[0], evalZ[-1]], 'b-', linewidth=2)
ax3.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)

ax3.set_xlim([0, xEnd])
ax3.set_ylim([np.min(evalZ), np.max(evalZ)])
ax3.xaxis.set_ticklabels([])
ax3.yaxis.set_ticklabels([])

imgname = imgdir + "volume-fraction-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##
 # Autocorrelate -- time slices
 ##
autoFig = plt.figure(figsize=(4,6))
gs = matplotlib.gridspec.GridSpec(3,1)

# autocorr of one location in time
bx1 = autoFig.add_subplot(gs[0,0])
vfAutoCorrSlice = AutoCorrelationFFT(vFrac[zz,:])
bx1.plot(time, vfAutoCorrSlice, 'k-', linewidth=2)

bx1.set_xlim([0, xEnd])
bx1.set_ylim([-.5, 1])
bx1.set_ylabel(r'$\langle \phi(t) \phi(t + \tau) \rangle$', fontsize=14)
bx1.xaxis.set_ticklabels([])
bx1.set_yticks([-0.5, 0 , 0.5, 1])

bx2 = autoFig.add_subplot(gs[1:3,0])

# autocorr of all locations in time
vfAutoCorr = np.zeros((nz, nt))
for zi,_ in enumerate(evalZ):
  # length of result is ceil(length(time)/2)
  vfAutoCorr[zi,:] = AutoCorrelationFFT(vFrac[zi,:])

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]])

bx2.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)

bx2.set_xlim([0, xEnd])
bx2.set_ylim([evalZ[0], evalZ[-1]])
bx2.set_xlabel(r"$\tau \ [s]$", fontsize=14)
bx2.set_ylabel(r'$z\ [mm]$', fontsize=14)

imgname = imgdir + "autocorrelate-time-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##
 # AutoCorrelate -- space slices
 ##
autoFig = plt.figure(figsize=(6,4))
gs = matplotlib.gridspec.GridSpec(1,3)

# Autocorrelate of one slice at constant time
cx3 = autoFig.add_subplot(gs[0,0])
vfAutoCorrSpaceSlice = AutoCorrelationFFT(vFrac[:,tt])

cx3.plot(vfAutoCorrSpaceSlice, dz, 'b-', linewidth=2)

cx3.set_xlim([-.5, 1])
cx3.set_ylim([dz[0], dz[-1]])
cx3.set_xlabel(r'$\langle \phi(z) \phi(z + \Delta z) \rangle$', fontsize=14)
cx3.set_ylabel(r'$\Delta z\ [mm]$', fontsize=14)
cx3.set_xticks([-0.5, 0 , 0.5, 1])

# autocorr space, all time
cx4 = autoFig.add_subplot(gs[0,1:3])
vfAutoCorr = np.zeros((nz, nt))
for ti,tval in enumerate(time):

  # length of result is ceil(length(time)/2)
  vfAutoCorr[:,ti] = AutoCorrelationFFT(vFrac[:,ti])

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])

cx4.plot([time[tt], time[tt]], [dz[0], dz[-1]], 'b-', linewidth=2)

cx4.set_xlim([0, xEnd])
cx4.set_ylim([dz[0], dz[-1]])
cx4.set_xlabel(r"$t\ [s]$", fontsize=14)
cx4.yaxis.set_ticklabels([])

imgname = imgdir + "autocorrelate-space-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##
 # cross-correlate
 ##
zz = np.floor(np.size(evalZ)/4) - 20
zz2 = zz + 150
crossCorrFig = plt.figure(figsize=(6,4.))
crossCorrFig.subplots_adjust(bottom=-0.75)

# location of slices to be xcorr'd
gs = matplotlib.gridspec.GridSpec(3,1)
dx1 = crossCorrFig.add_subplot(gs[0,:])

dx1.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
      extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
      vmin=0.3, vmax=0.4)

dx1.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)
dx1.plot([time[0], time[-1]], [evalZ[zz2], evalZ[zz2]], 'k--', linewidth=2)

dx1.set_xlim([0, xEnd])
dx1.set_ylim([np.min(evalZ), np.max(evalZ)])
dx1.set_xlabel(r'$t\ [s]$', fontsize=16)
dx1.set_ylabel(r'$z\ [mm]$',rotation=0, fontsize=16)
dx1.yaxis.set_label_coords(labelx, 0.5)
dx1.xaxis.set_ticklabels([])

# xcorr of two slices
dx2 = crossCorrFig.add_subplot(gs[1,:])
vfCrossCorrSlice = CrossCorrelationFFT(vFrac[zz,:], vFrac[zz2,:])

dx2.plot(time, vfCrossCorrSlice, 'k-', linewidth=2)

dx2.set_xlim([0, xEnd])
dx2.set_ylim([-0.20, 0.4])
dx2.set_ylabel(r'$\langle \phi(x, t) \phi(x + \Delta x, t + \tau)\rangle$', fontsize=16)
dx2.yaxis.set_label_coords(labelx, 0.5)
dx2.set_yticks([-0.20, 0.0, 0.2, 0.4])
dx2.xaxis.set_ticklabels([])

# xcorr of all slices
dx3 = crossCorrFig.add_subplot(gs[2,:])

vfCrossCorr = np.zeros((nz,nt))
for zi, zval in enumerate(evalZ):
  # correctly loop through domain
  if zz + zi >= nz:
    zInd = zz + zi - nz
  else:
    zInd = zz + zi
  vfCrossCorr[zi,:] = CrossCorrelationFFT(vFrac[zz,:], vFrac[zInd,:])

vfCrossCorr /= vfCrossCorr[0,0]

plt.imshow(vfCrossCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]])

dx3.plot([time[0], time[-1]], [dz[zz2] - dz[zz],dz[zz2] - dz[zz]], 'k-', 
  linewidth=2)

dx3.set_xlabel(r'$\tau\  [s]$', fontsize=16)
dx3.set_ylabel(r'$\Delta z\ [mm]$',rotation=0, fontsize=16)
dx3.yaxis.set_label_coords(labelx, 0.5)

dx3.set_xlim([0, xEnd])
dx3.set_ylim([dz[0], dz[-1]])

imgname = imgdir + "crosscorrelate-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
