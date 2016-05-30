#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Plotting Utility ---- "
print "      Radon Xform Plot"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)

nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])

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

# Read wavespeed data
wavedir = root + "simdata/waveData"
wavedata = np.genfromtxt(wavedir, skip_header=1)

for nn in np.arange(0,np.size(wavedata,0)):
  currN = wavedata[nn,0]
  currRho = wavedata[nn,1]
  if (nparts == currN) & (rho == currRho):
    wavespeed = wavedata[nn,2]
print "\n      Wavespeed = %.2f" % wavespeed

# Check for something gone wrong
if wavespeed == -1:
  print "Something is not quite right..."
  sys.exit()

# Pull output data -- each column is a different time
savefile = datadir + "radon-xform"
radonXForm = np.genfromtxt(savefile)

# Find limits
totZ = np.abs(evalZ[0]) + np.abs(evalZ[-1])
totT = np.abs(time[0]) + np.abs(time[-1])
cLength = np.sqrt(totZ**2 + totT**2)

# Plot
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
maxval = np.max(np.abs(radonXForm))
ax1.imshow(radonXForm, extent=(-90, 90, -cLength/2., cLength/2.), aspect='auto',
  vmin=-maxval, vmax=maxval)

ax1.set_xlabel(r'Wavespeed $[mm/s]$')
ax1.set_ylabel(r'Distance from center $[\sqrt{mm^2 + s^2}]$')
ax1.set_xlim([-90,90])
ax1.set_ylim([-cLength/2., cLength/2.])

# Stupid ticks
cTicks = np.linspace(-250, 250, 11).astype(np.int)
degTicks = np.rad2deg(np.arctan(cTicks/dzdt))
ax1.set_xticks(degTicks)
cTickLabels = cTicks.astype(np.str)
cTickLabels = ['250','','','100',
               '50','0','50',
               '100','','','250']
               # grrr
print cTickLabels
ax1.set_xticklabels(cTickLabels)
#ax1.set_xticks([-72, -36, 0, 36, 72])
#currTicks = ax1.xaxis.get_ticklocs()
#scaledTicks = np.tan(-np.deg2rad(currTicks)) * dzdt
#scaledTicks = np.round(scaledTicks*10.)/10
#ax1.set_xticklabels(scaledTicks)

# plot wavespeed from correlations
cloc = np.rad2deg(np.arctan(-wavespeed/dzdt))
ax1.plot([cloc, cloc], [-cLength/2., cLength/2.], 'k--')

imgname = imgdir + "radon-xform"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
