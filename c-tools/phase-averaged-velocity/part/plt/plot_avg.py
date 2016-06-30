#!/usr/bin/env python2
import sys, os, csv
import matplotlib
#matplotlib.use('Agg')      ## marcc savefig
matplotlib.use(u'Qt4Agg')  ## marcc plt.show
#matplotlib.use(u'GTKAgg')  ## lucan plt.show
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.fftpack as scifft

os.system('clear')


#print ""
#print " ---- Phase-Averaged Particle Velocity Plotting Utility ---- "
#print ""

# Parse command line args and set up directory structure
if len(sys.argv) > 2:
  simdir = sys.argv[1]
  tstart = float(sys.argv[2])
else:
  simdir = raw_input("      Simulation directory: ")
  tstart = float(raw_input("      Starting time [ms]: "))
  # TODO if tstart is -1 or empty, choose statsimtime

if not simdir.endswith('/'):
  simdir = simdir + '/'

home = os.path.expanduser("~")
root = home + "/scratch/triply_per/" + simdir

fdatadir = root + simdir + "flow_vel/data/"
pdatadir = root + simdir + "part_vel/data/"

# Find nparts and rho
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])

# Print simulation data
#print "      Sim root directory set to: " + root
#print "      Sim directory set to: " + simdir

# Set up output file paths
fdatadir = root + "flow_vel/data/"
fdataFile = fdatadir + "phaseAveragedFlowVel"
pdatadir = root + "part_vel/data/"
pdataFile = pdatadir + "particleAvgVel"
termFile = home + "/scratch/triply_per/simdata/singlePartSedi"

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(fdatadir):
  print "      " + fdatadir + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(pdatadir):
  print "      " + pdatadir + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Pull fluid velocity and time
ftime = np.genfromtxt(fdataFile, skip_header=1, usecols=0)
uf = np.genfromtxt(fdataFile, skip_header=1, usecols=1)
vf = np.genfromtxt(fdataFile, skip_header=1, usecols=2)
wf = np.genfromtxt(fdataFile, skip_header=1, usecols=3)

# Pull particle velocity and time
ptime = np.genfromtxt(pdataFile, skip_header=1, usecols=0)
up = np.genfromtxt(pdataFile, skip_header=1, usecols=1)
vp = np.genfromtxt(pdataFile, skip_header=1, usecols=2)
wp = np.genfromtxt(pdataFile, skip_header=1, usecols=3)

# Standard deviation of particle velocity (sdev AT a time step)
upSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=4)
vpSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=5)
wpSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=6)

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
commonMaxTime = np.min([ptime[-1], ftime[-1]])
time = np.arange(0., commonMaxTime + 0.00001, 2.)
#tIn = float(raw_input('Choose a time where you want to take a mean from: '))
tIn = tstart
ind = np.argwhere(time >= tIn)
#print "Time averaged over %f to %f" % (time[ind[0]], time[-1])


up = np.interp(time, ptime, up)
vp = np.interp(time, ptime, vp)
wp = np.interp(time, ptime, wp)

upSdev= np.interp(time, ptime, upSdev)
vpSdev= np.interp(time, ptime, vpSdev)
wpSdev= np.interp(time, ptime, wpSdev)

uf = np.interp(time, ftime, uf)
vf = np.interp(time, ftime, vf)
wf = np.interp(time, ftime, wf)

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

# particle / velocity ratio
print ""
print "max(wp/wf) = %.5f" % np.max(wp[ind]/wf[ind])
print ""


# Plot up,uf separately
phaseVel = plt.figure(figsize=(12,8))
phaseVel.suptitle('Particle Averaged Relative Vel', fontsize=20)

# u
uAx = phaseVel.add_subplot(311)
uAx.plot(time, urel, 'k-', linewidth=2)
#uAx.plot(time, up, 'k-', linewidth=2)
#uAx.plot(time, up + upSdev , 'k--', linewidth=2)
#uAx.plot(time, up - upSdev , 'k--', linewidth=2)
#uAx.plot(time, uf, 'b-', linewidth=2)
uAx.set_ylabel('u [m/s]')

# v
vAx = phaseVel.add_subplot(312)
vAx.plot(time, vrel, 'k-', linewidth=2)
#vAx.plot(time, vp, 'k-', linewidth=2)
#vAx.plot(time, vp + vpSdev , 'k--', linewidth=2)
#vAx.plot(time, vp - vpSdev , 'k--', linewidth=2)
#vAx.plot(time, vf, 'b-', linewidth=2)
vAx.set_ylabel('v [m/s]')

# w
wAx = phaseVel.add_subplot(313)
wAx.plot(time, wrel, 'k-', linewidth=2)
#wAx.plot(time, wp, 'k-', linewidth=2)
#wAx.plot([time[0], time[-1]], [0,0], 'k-', linewidth=2)
#wAx.plot(time, wp + wpSdev , 'k--', linewidth=2)
#wAx.plot(time, wp - wpSdev , 'k--', linewidth=2)
#wAx.plot(time, wf, 'b-', linewidth=2)
wAx.set_ylabel('w [m/s]')
wAx.set_xlabel('Time [s]')

#plt.show()

# Plot sdev/urel
fig2 = plt.figure(figsize=(12,8))
fig2.suptitle('Relative Velocity Mag', fontsize=20)
ax1 = fig2.add_subplot(111)
#ax1.plot(time[1:], up[1:]/uf[1:], 'k')
#ax1.plot(time[1:], vp[1:]/vf[1:], 'b')
#ax1.plot(time[1:], wp[1:]/wf[1:], 'g')
ax1.set_ylim([-0.01, 0.01])

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

