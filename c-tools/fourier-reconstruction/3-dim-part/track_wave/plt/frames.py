#!/usr/bin/env python2
from setup import *
from functions import *
os.system('clear')

print ""
print " ---- 3D Fourier Reconstruction Plotting Utility ---- "
print "                Moving Average"
print ""

# Simulation size
domX = 42
domY = 42
domZ = 126

# Setup simulation parameters, directory structures, and get data
(partR, simdir) = simParams(sys)
(root, simdir, datadir, imgdir) = directoryStructure(simdir)
(time, nt, evalZ, nz) = initData(datadir)
printSimulationData(partR, root, simdir, datadir)

# Get mean volume fraction
nparts = int(simdir.partition('/')[0])
avgPhi = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)

# sort the moving_files
moving_files = sorted_nicely(glob.glob(datadir + "vfrac-moving-avg-*"))
stationary_files = sorted_nicely(glob.glob(datadir + "vfrac-stationary-avg-*"))
nFiles = len(moving_files)

print "      Found %d moving files" % nFiles

# Set up figure
fig = plt.figure(figsize=(5,5))
gs = matplotlib.gridspec.GridSpec(3,2)

moving_avg = np.genfromtxt(moving_files[0])
stat_avg = np.genfromtxt(stationary_files[0])

ax1 = fig.add_subplot(gs[0,0])
cInit = ax1.imshow(moving_avg, origin="lower", aspect="auto", interpolation="none",
  extent=[-0.5*domX, 0.5*domX, -0.5*domY, 0.5*domY], cmap ='seismic',
  vmin=0, vmax=2.*avgPhi)
ax1.set_xlabel("x")
ax1.set_xlim([-0.5*domX, 0.5*domX])
ax1.set_ylabel("y")
ax1.set_ylim([-0.5*domY, 0.5*domY])

plt.colorbar(cInit)

plt.cla()
  
ax2 = fig.add_subplot(gs[0,1])
cInit = ax2.imshow(stat_avg, origin="lower", aspect="auto", interpolation="none",
  extent=[-0.5*domX, 0.5*domX, -0.5*domY, 0.5*domY], cmap ='seismic',
  vmin=0, vmax=2.*avgPhi)
ax2.set_xlabel("x")
ax2.set_xlim([-0.5*domX, 0.5*domX])
ax2.set_ylabel("y")
ax2.set_ylim([-0.5*domY, 0.5*domY])

plt.colorbar(cInit)

plt.cla()

ax3 = fig.add_subplot(gs[1:3,0:2])
vfrac = np.genfromtxt(root + simdir + "../f-rec-1D/data/volume-fraction").T
time_vf = np.genfromtxt(root + simdir + "../f-rec-1D/data/info", skip_footer=1)[1:] / 1000
evalZ_vf = np.genfromtxt(root + simdir + "../f-rec-1D/data/info", skip_header=1)[1:]

minVal = np.floor(100*np.amin(vfrac))/100
maxVal = np.ceil(100*np.amax(vfrac))/100
upperDiff = maxVal - avgPhi
lowerDiff = avgPhi - minVal
maxDiff = np.max([upperDiff, lowerDiff])

ax3.imshow(vfrac, origin="lower", aspect="auto", interpolation="none",
  extent=[time_vf[0], time_vf[-1], evalZ_vf[0], evalZ_vf[-1]],
  vmin=avgPhi - maxDiff, vmax=avgPhi + maxDiff,
  cmap='seismic')

yerr = 1.5*partR
ax3.errorbar(0.4, evalZ[0], fmt='ko', yerr=yerr)

ax3.set_xlabel('time')
ax3.set_xlim([.400, .600])
ax3.set_ylabel('z')
ax3.set_ylim([-5, 25])

#imgname = imgdir + "test"
#plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#sys.exit()

plt.sca(ax1)

# Loop over moving_files and save as png
for ff, fname in enumerate(moving_files):
  moving_avg = np.genfromtxt(fname)
  stat_avg = np.genfromtxt(stationary_files[ff])

  ax1.imshow(moving_avg, origin="lower", aspect="auto", interpolation="none",
    extent=[-0.5*domX, 0.5*domX, -0.5*domY, 0.5*domY], cmap ='seismic',
    vmin=0, vmax=2.*avgPhi)
  ax1.set_title("time = %.3f, z = %.2f" % (time[ff], evalZ[ff]))

  plt.sca(ax2)
  ax2.imshow(stat_avg, origin="lower", aspect="auto", interpolation="none",
    extent=[-0.5*domX, 0.5*domX, -0.5*domY, 0.5*domY], cmap ='seismic',
    vmin=0, vmax=2.*avgPhi)
  ax2.set_title("time = %.3f, z = %.2f" % (time[ff], evalZ[0]))

  plt.sca(ax3)
  ax3.errorbar(time[ff], evalZ[ff], fmt='ko', yerr=yerr)
  ax3.errorbar(time[ff], evalZ[0], fmt='ko', yerr=yerr)

  #imgname = imgdir + fname.split('/')[-1]
  imgname = imgdir + "vfrac-avg-%05d" % ff
  plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

  plt.sca(ax2)
  plt.cla()
  plt.sca(ax1)
  plt.cla()
