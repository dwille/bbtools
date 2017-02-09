#!/usr/bin/env python2
from setup import *

os.system('clear')
print ""
print " ---- Phase-Averaged Particle Velocity Plotting Utility ---- "
print ""

# Parse command line args and set up directory structure
(partR, simdir, tstart) = simParams(sys)
(root, simdir, pdatadir, fdatadir, imgdir) = directoryStructure(simdir)
printSimulationData(partR, root, simdir, pdatadir, fdatadir)

# Find nparts and rho
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][3:6])

# Set up output file paths
fdataFile = fdatadir + "phaseAveragedFlowVel"
pdataFile = pdatadir + "particleAvgVel"
termFile = root + "simdata/singlePartSedi"

# Check if terminal vel file exists
if not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Pull fluid velocity and time AFTER steady state time
ftime = np.genfromtxt(fdataFile, skip_header=1, usecols=0)
tf_start_ind = np.squeeze(np.argwhere(ftime >= tstart)[0])

uf = np.genfromtxt(fdataFile, skip_header=1, usecols=1)[tf_start_ind:]
vf = np.genfromtxt(fdataFile, skip_header=1, usecols=2)[tf_start_ind:]
wf = np.genfromtxt(fdataFile, skip_header=1, usecols=3)[tf_start_ind:]

# Pull particle velocity and time AFTER steady state time
ptime = np.genfromtxt(pdataFile, skip_header=1, usecols=0)
tp_start_ind = np.squeeze(np.argwhere(ptime >= tstart)[0])

up = np.genfromtxt(pdataFile, skip_header=1, usecols=1)[tp_start_ind:]
vp = np.genfromtxt(pdataFile, skip_header=1, usecols=2)[tp_start_ind:]
wp = np.genfromtxt(pdataFile, skip_header=1, usecols=3)[tp_start_ind:]

# Standard deviation of particle velocity (sdev AT a time step)
upSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=4)[tp_start_ind:]
vpSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=5)[tp_start_ind:]
wpSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=6)[tp_start_ind:]

# Pull terminal velocity, find correct for the simulation
termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)
if rho == 2.0:
  termVel = termVel[0]
elif rho == 3.3:
  termVel = termVel[1]
elif rho == 4.0:
  termVel = termVel[2]
elif rho == 5.0:
  termVel = termVel[3]

# Interpolate fluid/particle times together
commonMinTime = np.max([ptime[0], ftime[0]);
commonMaxTime = np.min([ptime[-1], ftime[-1]])
dt = np.mean(np.diff(ptime))
time = np.arange(commonMinTime - 1e-5, commonMaxTime + 1e-5, dt)

up = np.interp(time, ptime, up)
vp = np.interp(time, ptime, vp)
wp = np.interp(time, ptime, wp)

uf = np.interp(time, ftime, uf)
vf = np.interp(time, ftime, vf)
wf = np.interp(time, ftime, wf)

upSdev= np.interp(time, ptime, upSdev)
vpSdev= np.interp(time, ptime, vpSdev)
wpSdev= np.interp(time, ptime, wpSdev)

# Relative velocity
urel = uf - up
vrel = vf - vp
wrel = wf - wp

# Normalize by relative velocity
ufluct_rel = np.zeros(np.size(time))
vfluct_rel = np.zeros(np.size(time))
wfluct_rel = np.zeros(np.size(time))
ufluct_rel[1:] = upSdev[1:] / wrel[1:]
vfluct_rel[1:] = vpSdev[1:] / wrel[1:]
wfluct_rel[1:] = wpSdev[1:] / wrel[1:]

# Normalize by terminal velocity
ufluct_term = upSdev / termVel
vfluct_term = vpSdev / termVel
wfluct_term = wpSdev / termVel

# For averaging velocities over time
#tIn = float(raw_input('Choose a time where you want to take a mean from: '))
tIn = tstart
ind = np.argwhere(time >= tIn)
#print "Time averaged over %f to %f" % (time[ind[0]], time[-1])

# particle / fluid velocity ratio
# print ""
# print "max(wp/wf) = %.5f" % np.max(wp[ind]/wf[ind])
# print ""

# Plot urel
uFig = plt.figure()

# urel
ax1 = uFig.add_subplot(111)
ax1.plot(time, urel, 'k-', linewidth=2)
ax1.plot(time, vrel, 'b-', linewidth=2)
ax1.plot(time, wrel, 'r-', linewidth=2)

ax1.set_ylabel(r'$\langle u_i^f \rangle - \langle u_i^p \rangle\ (m/s)$')
ax1.set_xlabel(r'$t\ (s)$')

ax1.legend([r"$u$", r"$v$", r"$w$"])


#plt.show()
imgname = imgdir + "relative_vel"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

#uAx.plot(time, up, 'k-', linewidth=2)
#uAx.plot(time, up + upSdev , 'k--', linewidth=2)
#uAx.plot(time, up - upSdev , 'k--', linewidth=2)
#uAx.plot(time, uf, 'b-', linewidth=2)
#vAx.plot(time, vp, 'k-', linewidth=2)
#vAx.plot(time, vp + vpSdev , 'k--', linewidth=2)
#vAx.plot(time, vp - vpSdev , 'k--', linewidth=2)
#vAx.plot(time, vf, 'b-', linewidth=2)
#wAx.plot(time, wp, 'k-', linewidth=2)
#wAx.plot([time[0], time[-1]], [0,0], 'k-', linewidth=2)
#wAx.plot(time, wp + wpSdev , 'k--', linewidth=2)
#wAx.plot(time, wp - wpSdev , 'k--', linewidth=2)
#wAx.plot(time, wf, 'b-', linewidth=2)

# # Plot sdev/urel
# fig2 = plt.figure(figsize=(12,8))
# fig2.suptitle('Relative Velocity Mag', fontsize=20)
# ax1 = fig2.add_subplot(111)
# #ax1.plot(time[1:], up[1:]/uf[1:], 'k')
# #ax1.plot(time[1:], vp[1:]/vf[1:], 'b')
# #ax1.plot(time[1:], wp[1:]/wf[1:], 'g')
# ax1.set_ylim([-0.01, 0.01])

#plt.show()

## Normalized fluctuations -- u_sdev / u_rel
#print "Steady state value of u/urel = %.5f %.5f %.5f" % \
#  (np.mean(ufluct_rel[ind]), np.mean(vfluct_rel[ind]), 
#  np.mean(wfluct_rel[ind]))

## Normalized fluctuations -- u_sdev / u_term
#print "Steady state value of u/u_term = %.5f %.5f %.5f" % \
#  (np.mean(ufluct_term[ind]), np.mean(vfluct_term[ind]), 
#  np.mean(wfluct_term[ind]))

## Raw fluctuations
## steady state u',v',w',urel,vrel,rel
#print "Mean of up_sdev, u_rel after given time:"
#print "upSdev, vpSdev, wpSdev, urel, vrel, wrel"
#print "%.5f %.5f %.5f %.5f %.5f %.5f" % \
#  (np.mean(upSdev[ind]), np.mean(vpSdev[ind]), np.mean(wpSdev[ind]), 
#   np.mean(urel[ind]), np.mean(vrel[ind]), np.mean(wrel[ind]))

## Mean of uf - up
#print "Mean of uf - up after given time:"
#print "%.5f, %.5f, %.5f" % \
#  (np.mean(urel[ind]), np.mean(vrel[ind]), np.mean(wrel[ind]))

## Mean of uf
#print "Mean of uf after given time:"
#print "%.5f, %.5f, %.5f" % \
#  (np.mean(uf[ind]), np.mean(vf[ind]), np.mean(wf[ind]))

