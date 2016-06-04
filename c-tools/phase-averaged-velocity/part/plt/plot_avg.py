#!/usr/bin/env python2

import matplotlib
import matplotlib.rcsetup as rcsetup
matplotlib.use(u'Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
os.system('clear')

## GET INFO
print " ---- Phase-Averaged Particle Velocity Plotting Utility ---- "
print ""

# DEVEL
#root = "/home/dwille/bbtools/c-tools/phase-averaged-velocity/"
#simdir = "sim/"
#datadir = root + simdir + "data-reconstruct/"

# MARCC
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
simdir = raw_input("      Simulation directory: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'
fdatadir = root + simdir + "flow_vel/data/"
pdatadir = root + simdir + "part_vel/data/"

print "      Sim root directory set to: " + root
print "      Using location: " + fdatadir
print "      Using location: " + pdatadir

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(fdatadir):
  print "      " + fdatadir + " does not exist. Exiting..."
  print ""
  sys.exit()
  
if not os.path.exists(pdatadir):
  print "      " + pdatadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + simdir + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Set up output file paths
fdataFile = fdatadir + "phaseAveragedFlowVel"
pdataFile = pdatadir + "particleAvgVel"

ftime = np.genfromtxt(fdataFile, skip_header=1, usecols=0)
uf = np.genfromtxt(fdataFile, skip_header=1, usecols=1)
vf = np.genfromtxt(fdataFile, skip_header=1, usecols=2)
wf = np.genfromtxt(fdataFile, skip_header=1, usecols=3)

ptime = np.genfromtxt(pdataFile, skip_header=1, usecols=0)
up = np.genfromtxt(pdataFile, skip_header=1, usecols=1)
vp = np.genfromtxt(pdataFile, skip_header=1, usecols=2)
wp = np.genfromtxt(pdataFile, skip_header=1, usecols=3)
upSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=4)
vpSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=5)
wpSdev = np.genfromtxt(pdataFile, skip_header=1, usecols=6)

# Interpolate fluid/particle times together
commonMaxTime = np.min([ptime[-1], ftime[-1]])
time = np.arange(0., commonMaxTime + 0.00001, 2.)

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

# Magnitude of fluctuations
#ufluct = upSdev / urel
#vfluct = vpSdev / vrel
wfluct = wpSdev / wrel

# Plot up,uf separately
phaseVel = plt.figure(figsize=(12,8))
phaseVel.suptitle('Particle Averaged Vel', fontsize=20)

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

plt.show()

# Plot sdev/urel
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
#ax1.plot(time, ufluct, 'k')
#ax1.plot(time, vfluct, 'b')
ax1.plot(time, wfluct, 'g')

plt.show()

tIn = float(raw_input('Choose a time where you want to take a mean from: '))
ind = np.argwhere(time >= tIn)
print "Steady state value of ufluct = %.5f" % np.mean(ufluct[ind])
print "Steady state value of vfluct = %.5f" % np.mean(vfluct[ind])
print "Steady state value of wfluct = %.5f" % np.mean(wfluct[ind])
#print "Steady state value of up = %.5f" % np.mean(up[ind])
#print "Steady state value of vp = %.5f" % np.mean(vp[ind])
#print "Steady state value of wp = %.5f" % np.mean(wp[ind])
#print "Steady state sdev of up = %.5f" % np.mean(upSdev[ind])
#print "Steady state sdev of vp = %.5f" % np.mean(vpSdev[ind])
#print "Steady state sdev of wp = %.5f" % np.mean(wpSdev[ind])
