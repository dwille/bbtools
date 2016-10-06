#!/usr/bin/env python2
from setup import *
os.system('clear')

from matplotlib.ticker import MultipleLocator

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "                Slices Plotting"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])
nu = .01715   ## mm^2/ms XXX

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

imgdir = imgdir + "slices/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)
dz = evalZ - evalZ[0]

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# get avg volume fraction
domX = 42   # XXX
domY = 42   # XXX
domZ = 126  # XXX
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

##
 # Read Volume Fraction
 ##
vFracFile = datadir + "volume-fraction"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]

# Slice location
tt = 250
zz = np.floor(np.size(evalZ)/2)
labelx = -0.15

#tau = (2.*partR) / phaseVel       ## mm / (mm/s) = s
# Non-dimensionalize
evalZ /= (2.*partR)               ## mm / mm
dz /= (2.*partR)               ## mm / mm
tau = (2.*partR)*(2.*partR)/nu    ## mm^2/(mm^2/ms) = ms
tau /= 1000                       ## ms -> s
time /= tau                       ## s / s
xEnd = time[-1]
xEnd = 100.*np.floor(xEnd / 100.)

# minmax of vfrac
minVal = np.floor(100*np.amin(vFrac))/100
maxVal = np.ceil(100*np.amax(vFrac))/100
upperDiff = maxVal - avgVolumeFraction
lowerDiff = avgVolumeFraction - minVal
maxDiff = np.max([upperDiff, lowerDiff])
maxDiff = np.floor(maxDiff * 100.) / 100.

# ##
#  # Plot one volume fraction slice
#  ##
# vFracFig = plt.figure(figsize=(6,4))
# gs = matplotlib.gridspec.GridSpec(2,3)
# 
# 
# # fixed time slice
# ax1 = vFracFig.add_subplot(gs[0,0])
# ax1.plot(vFrac[:,tt], evalZ, 'k-', linewidth=2)
# 
# ax1.set_xlabel(r'$\phi$', fontsize=14)
# ax1.xaxis.set_major_locator(MultipleLocator(0.05))
# ax1.xaxis.set_minor_locator(MultipleLocator(0.025))
# ax1.set_xlim([0.3, 0.4])
# 
# ax1.set_ylabel(r'$\frac{z}{2a}$', fontsize=14, rotation=0)
# ax1.yaxis.set_major_locator(MultipleLocator(5))
# ax1.yaxis.set_minor_locator(MultipleLocator(2.5))
# ax1.set_ylim([np.min(evalZ), np.max(evalZ)])
# 
# #ax1.annotate(r"$(a)$",xy=get_axis_limits(ax1))
# 
# # fixed location slice
# ax2 = vFracFig.add_subplot(gs[1,1:3])
# plt.plot(time, vFrac[zz,:], 'k', linewidth=2)
# 
# plt.xlim([0, xEnd])
# plt.xticks(np.floor(np.arange(0, xEnd+0.01, 100)))
# 
# ax2.set_xlabel(r"$\frac{\langle w_f \rangle t}{2a}$", fontsize=14)
# ax2.xaxis.set_major_locator(MultipleLocator(100))
# ax2.xaxis.set_minor_locator(MultipleLocator(50))
# ax2.set_xlim([0, xEnd])
# 
# ax2.set_ylabel(r"$\phi$", fontsize=14)
# ax2.yaxis.set_major_locator(MultipleLocator(0.05))
# ax2.yaxis.set_minor_locator(MultipleLocator(0.025))
# ax2.set_ylim([0.3, 0.4])
# 
# #ax2.annotate(r"$(c)$",xy=get_axis_limits(ax2))
# 
# # contour plot
# ax3 = vFracFig.add_subplot(gs[0,1:3])
# 
# im = ax3.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
#       extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
#       vmin=0.3, vmax=0.4, cmap='coolwarm')
# 
# ax3.plot([time[tt], time[tt]], [evalZ[0], evalZ[-1]], 'k-', linewidth=2)
# ax3.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)
# #ax3.plot([time[0], time[-1]], [evalZ[zz+50], evalZ[zz+50]], 'k--', linewidth=2)
# 
# ax3.xaxis.set_major_locator(MultipleLocator(100))
# ax3.xaxis.set_minor_locator(MultipleLocator(50.))
# ax3.set_xlim([0, xEnd])
# 
# ax3.yaxis.set_major_locator(MultipleLocator(5))
# ax3.yaxis.set_minor_locator(MultipleLocator(2.5))
# ax3.set_ylim([np.min(evalZ), np.max(evalZ)])
# 
# ax3.xaxis.set_ticklabels([])
# ax3.yaxis.set_ticklabels([])
# 
# #box_props = dict(boxstyle="circle,pad=0.3", fc="w", ec="w", lw=0)
# #ax3.annotate(r"$\ (b)$",xy=get_axis_limits(ax3), bbox=box_props)
# 
# imgname = imgdir + "volume-fraction-slices"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
# plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##
 # cross-correlate
 ##
zz = np.floor(np.size(evalZ)/4) + 10
zz2 = zz + 150

# location of slices to be xcorr'd
vFracFig = plt.figure(figsize=(3.25,1.625))
plt.imshow(vFrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=avgVolumeFraction - maxDiff, vmax=avgVolumeFraction + maxDiff,
  cmap='coolwarm')

plt.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)
plt.plot([time[0], time[-1]], [evalZ[zz2], evalZ[zz2]], 'k--', linewidth=2)

# XXX
#plt.annotate(s='', xy=(0,-7.5), xytext=(0.3,-7.5), arrowprops=dict(arrowstyle='-'))
#plt.annotate(s='', xy=(0.3,-7.4), xytext=(0.3,evalZ[zz2]), arrowprops=dict(arrowstyle='-'))


#cbar = plt.colorbar()
plt.xlabel(r"$\nu t/(2a)^2$")
plt.ylabel(r'$z/2a$', rotation=90)
plt.gca().yaxis.set_label_coords(-0.14, 0.45)
plt.xlim([0, 6])  # XXX
plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(5))
plt.gca().yaxis.set_minor_locator(MultipleLocator(2.5))

imgname = imgdir + "volume-fraction-xcorr-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# xcorr of two slices
fig2 = plt.figure(figsize=(3.25,1.625))
vfAutoCorr = CrossCorrelationFFT(vFrac[zz,:], vFrac[zz,:])
vfCrossCorrSlice = CrossCorrelationFFT(vFrac[zz,:], vFrac[zz2,:])
vfCrossCorrSlice /= vfAutoCorr[0]

plt.plot(time, vfCrossCorrSlice, 'k-', linewidth=2)

plt.xlabel(r"$\nu t/(2a)^2$")
#plt.ylabel(r'$\langle \phi(z_*, t_*) \phi(z_* + \Delta z_*, t_* + \Delta t_*) \rangle$', rotation=90)
plt.ylabel(r'$(\phi \star \phi)(\Delta z_*, \Delta t_*)$', rotation=90)
plt.gca().yaxis.set_label_coords(-0.14, 0.45)

plt.xlim([0, 6])
plt.ylim([-1,1])
plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
plt.gca().yaxis.set_minor_locator(MultipleLocator(0.25))

imgname = imgdir + "volume-fraction-xcorr"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# xcorr of all slices
fig3 = plt.figure(figsize=(3.25,1.625))

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
  extent=[time[0], time[-1], dz[0], dz[-1]], 
  vmin=-1., vmax=1., cmap='coolwarm')

#plt.plot([time[0], time[-1]], [dz[zz2] - dz[zz],dz[zz2] - dz[zz]], 'k-', 
#  linewidth=2)

plt.xlabel(r'$\Delta t_*$')
plt.xlim([0, 6])
plt.ylabel(r'$\Delta z_*$', rotation=0)
plt.ylim([dz[0], dz[-1]])
plt.gca().yaxis.set_label_coords(-0.15, 0.45)

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(5))
plt.gca().yaxis.set_minor_locator(MultipleLocator(2.5))

imgname = imgdir + "crosscorrelate-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##
 # Plot one velocity slice
 ##
wpFile = datadir + "part-w"
wp = np.genfromtxt(wpFile).T[:,tsInd:]
wpFig = plt.figure(figsize=(6,4))
gs = matplotlib.gridspec.GridSpec(2,3)

# fixed time slice
ax1 = wpFig.add_subplot(gs[0,0])
ax1.plot(wp[:,tt], evalZ, 'k-', linewidth=2)

ax1.set_xlabel(r'$w_p$', fontsize=14)
ax1.xaxis.set_major_locator(MultipleLocator(0.025))
ax1.xaxis.set_minor_locator(MultipleLocator(0.0125))
ax1.set_xlim([-0.025, 0.025])

ax1.set_ylabel(r'$z\ [\mathrm{mm}]$', fontsize=14)
ax1.yaxis.set_major_locator(MultipleLocator(20))
ax1.yaxis.set_minor_locator(MultipleLocator(10))
ax1.set_ylim([np.min(evalZ), np.max(evalZ)])

ax1.annotate(r"$(a)$",xy=get_axis_limits(ax1))

# fixed location slice
ax2 = wpFig.add_subplot(gs[1,1:3])
plt.plot(time, wp[zz,:], 'k', linewidth=2)

xEnd = time[-1]

ax2.set_xlabel(r"$t\ [\mathrm{s}]$", fontsize=14)
ax2.xaxis.set_major_locator(MultipleLocator(1))
ax2.xaxis.set_minor_locator(MultipleLocator(0.5))
plt.xlim([0, xEnd])

ax2.set_ylabel(r"$w_p$", fontsize=14)
ax2.yaxis.set_major_locator(MultipleLocator(0.025))
ax2.yaxis.set_minor_locator(MultipleLocator(0.0125))
ax2.set_ylim([-0.025, 0.025])

ax2.annotate(r"$(c)$",xy=get_axis_limits(ax2))

# contour plot
ax3 = wpFig.add_subplot(gs[0,1:3])

im = ax3.imshow(wp, origin="lower", aspect="auto", interpolation="none",
      extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
      vmin=-0.025, vmax=0.025, cmap='seismic')

ax3.plot([time[tt], time[tt]], [evalZ[0], evalZ[-1]], 'k-', linewidth=2)
ax3.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)

ax3.set_xlim([0, xEnd])
ax3.xaxis.set_major_locator(MultipleLocator(1))
ax3.xaxis.set_minor_locator(MultipleLocator(0.5))

ax3.yaxis.set_major_locator(MultipleLocator(20))
ax3.yaxis.set_minor_locator(MultipleLocator(10))
ax3.set_ylim([np.min(evalZ), np.max(evalZ)])

ax3.xaxis.set_ticklabels([])
ax3.yaxis.set_ticklabels([])

box_props = dict(boxstyle="circle,pad=0.3", fc="w", ec="w", lw=0)
ax3.annotate(r"$\ (b)$",xy=get_axis_limits(ax3), bbox=box_props)

imgname = imgdir + "vertical-velocity-slices"
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
bx1.set_ylabel(r'$\langle \phi(%d,t) \phi(%d,t + \Delta t) \rangle$' % (evalZ[zz],evalZ[zz]),
  fontsize=14)
bx1.xaxis.set_ticklabels([])
bx1.set_yticks([-0.5, 0 , 0.5, 1])

bx2 = autoFig.add_subplot(gs[1:3,0])

# autocorr of all locations in time
vfAutoCorr = np.zeros((nz, nt))
for zi,_ in enumerate(evalZ):
  # length of result is ceil(length(time)/2)
  vfAutoCorr[zi,:] = AutoCorrelationFFT(vFrac[zi,:])

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]], 
  vmin=-1., vmax=1., cmap='seismic')

bx2.plot([time[0], time[-1]], [evalZ[zz], evalZ[zz]], 'k-', linewidth=2)

bx2.set_xlim([0, xEnd])
bx2.set_ylim([evalZ[0], evalZ[-1]])
bx2.set_xlabel(r"$\Delta t \ [\mathrm{s}]$", fontsize=14)
bx2.set_ylabel(r'$z\ [\mathrm{mm}]$', fontsize=14)

imgname = imgdir + "autocorrelate-time-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##
 # AutoCorrelate -- space slices
 ##
autoFig = plt.figure(figsize=(8,4))
gs = matplotlib.gridspec.GridSpec(1,3)

# Autocorrelate of one slice at constant time
cx3 = autoFig.add_subplot(gs[0,0])
vfAutoCorrSpaceSlice,_,_ = AutoCorrelationSpaceFFT(vFrac[:,tt])

cx3.plot(vfAutoCorrSpaceSlice, dz, 'b-', linewidth=2)

cx3.set_xlim([-.5, 1])
cx3.set_ylim([dz[0], dz[-1]])
cx3.set_xlabel(r'$\langle \phi(z) \phi(z + \Delta z) \rangle$', fontsize=14)
cx3.set_ylabel(r'$\Delta z\ [\mathrm{mm}]$', fontsize=14)
cx3.set_xticks([-0.5, 0 , 0.5, 1])

# autocorr space, all time
cx4 = autoFig.add_subplot(gs[0,1:3])
vfAutoCorr = np.zeros((nz, nt))
for ti,tval in enumerate(time):

  # length of result is ceil(length(time)/2)
  vfAutoCorr[:,ti],_,_ = AutoCorrelationSpaceFFT(vFrac[:,ti])

plt.imshow(vfAutoCorr, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], dz[0], dz[-1]], 
  vmin=-1., vmax=1., cmap='seismic')

cx4.plot([time[tt], time[tt]], [dz[0], dz[-1]], 'b-', linewidth=2)

cx4.set_xlim([0, xEnd])
cx4.set_ylim([dz[0], dz[-1]])
cx4.set_xlabel(r"$t\ [\mathrm{s}]$", fontsize=14)
cx4.yaxis.set_ticklabels([])

imgname = imgdir + "autocorrelate-space-slices"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
