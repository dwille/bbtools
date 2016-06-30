#!/usr/bin/env python2

from mpl_toolkits.mplot3d import axes3d
import matplotlib
matplotlib.use(u'Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import os,sys
os.system('clear')


print ""
print " ---- Phase-Averaged Particle Velocity Flucts Plotting Utility"
print ""

# Simulation Params
d = 4.2
nu = .01715
rhof = 8.75e-4
partR = 2.1
g = 0.00981
rhopf = np.array([2.0, 3.3, 4.0, 5.0])
rhop = np.array([0.001750, .00289, 0.0035, 0.004375])
phi = np.array([0.087, 0.175, 0.262, 0.349])

# MARCC
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
fluctFile = root + "velFlucts_raw"
termFile = root + "singlePartSedi"

print "      Root directory set to: " + root
print "      Using location: " + fluctFile
print "      Using location: " + termFile
print ""
print "      Using d = %.2f" % d
print "      Using nu = %.5f" % nu

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(fluctFile):
  print "      " + fluctFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

imgdir = imgdir + "velocity_flucts/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Read data
wflucts = np.genfromtxt(fluctFile, skip_header=1, usecols=4)
dwdphi = np.genfromtxt(fluctFile, skip_header=1, usecols=8)

termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)
wrel = np.genfromtxt(fluctFile, skip_header=1, usecols=7)

# normalize
#wflucts[0:4] /= termVel
#wflucts[4:8] /= termVel
#wflucts[8:12] /= termVel
#wflucts[12:16] /= termVel
#wflucts /= wrel

Re_t = d*termVel/nu
Ga = np.sqrt(np.abs(rhopf - 1.)*g*(2.*partR)**3.)/nu

## Plot as a function of phi
fig1 = plt.figure(figsize=(5,5))

ax = fig1.gca(projection="3d")

yax = Ga
# Surface
ax.plot(phi[0]*np.ones(4), yax, wflucts[0:4], 'b*', markersize=7)
ax.plot(phi[1]*np.ones(4), yax, wflucts[4:8], 'gs', markersize=7)
ax.plot(phi[2]*np.ones(4), yax, wflucts[8:12], 'ro', markersize=7)
ax.plot(phi[3]*np.ones(4), yax, wflucts[12:16], 'c^', markersize=7)

ax.set_xlabel(r"$\phi$")
ax.set_xlim([0.0,0.45])
ax.set_xticks([0, 0.15, 0.30, 0.45])

ax.set_ylabel(r"$Ga$")
ax.set_ylim([20,120])
ax.set_yticks([20,40,60,80,100,120])

ax.set_zlabel(r"$w^\prime$")
#ax.set_zlim([0.0,0.10])
#ax.set_zticks([0,.02,.04,.06,.08,.10])

xmin = np.min(ax.get_xlim())*np.ones(4)
ymax = np.max(ax.get_ylim())*np.ones(4)
# Projection -- w(Re_t)
ax.plot(phi[0]*np.ones(4), ymax, wflucts[0:4], 'b*', markersize=4, alpha=0.4)
ax.plot(phi[1]*np.ones(4), ymax, wflucts[4:8], 'gs', markersize=4, alpha=0.4)
ax.plot(phi[2]*np.ones(4), ymax, wflucts[8:12], 'ro', markersize=4, alpha=0.4)
ax.plot(phi[3]*np.ones(4), ymax, wflucts[12:16], 'c^', markersize=4, alpha=0.4)

#ax.plot(phi, np.ones(4), wflucts[0:16:4], 'b*', zdir='y')

# Projection -- w(phi)
ax.plot(xmin, yax, wflucts[0:4], 'b*', markersize=4, alpha=0.4)
ax.plot(xmin, yax, wflucts[4:8], 'gs', markersize=4, alpha=0.4)
ax.plot(xmin, yax, wflucts[8:12], 'ro', markersize=4, alpha=0.4)
ax.plot(xmin, yax, wflucts[12:16], 'c^', markersize=4, alpha=0.4)



plt.tight_layout()
imgname = imgdir + "velflucts_3d"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
