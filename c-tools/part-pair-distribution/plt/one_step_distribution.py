#!/usr/bin/env python2
import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt

from floating_polar import fractional_polar_axes

print ""
print " ---- Particle Pair Distribution Plotting Utility ---- "
print ""

## Nicely sorted function
def sorted_nicely( l ):
  """ Sorts the given iterable in a natural way

  Required arguments:
  l -- The iterable to be sroted.

  courtesy stackoverflow/2669059
  """
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
  return sorted(l, key = alphanum_key)

# SIMULATION PARAMETERS
partR = 2.1

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
datadir = root + "part-pair/data/"

# Find nparts and rho
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])

# Print simulation data
print "      Sim root directory set to: " + root
print "      Sim directory set to: " + simdir

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(datadir):
  print "      " + datadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Data Files
infoFile = datadir + "info"

# Create imgdir if necessary
imgdir = root + "part-pair/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=2)[1:]
nt = np.size(time)

evalR = np.genfromtxt(infoFile, skip_header=1,skip_footer=1)[1:]/(2.*partR)
nr = np.size(evalR)

evalTh = np.genfromtxt(infoFile, skip_header=2)[1:]
nth = np.size(evalTh)

data = np.genfromtxt(datadir + "part-pair-avg")

## Find output data -- each column is a different time
#files = sorted_nicely(glob.glob(datadir + "part-pair-*"))
#nFiles = len(files)
#
#print "      Found " + str(nFiles) + " files."
#
## Loop and pull
##data = np.zeros((nt, nr, nth));
#data = np.zeros((nr, nth));
#for ff, fname in enumerate(files):
#  #data[ff,:,:] = np.genfromtxt(fname)
#  data += np.genfromtxt(fname)/float(nFiles)

## Plot ###
#theta,rad = np.meshgrid(np.flipud(evalTh), evalR)
theta,rad = np.meshgrid(evalTh, evalR)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, polar = True)
#ax1.pcolormesh(theta,rad, data[0,:,:])
c = ax1.pcolormesh(theta,rad, data[:,:], norm=None)
cbaxes = fig1.add_axes([0.45, 0.55, 0.03, 0.4])
fig1.colorbar(c,ax=ax1, cax=cbaxes)
cbaxes.yaxis.set_ticks_position('left')

ax1.set_theta_zero_location("N")
ax1.set_theta_direction('clockwise')
ax1.set_xticks([0., np.pi/6., np.pi/3., np.pi/2.])
ax1.set_xticklabels([r"$0$", r"$\pi/6$", r"$\pi/3$", r"$\pi/2$"])
ax1.set_ylim([0, np.amax(evalR)])
ax1.set_rlabel_position(90.)
label_position = ax1.get_rlabel_position()
ax1.text(np.radians(label_position+10), ax1.get_rmax()/2., r'$r/d$',
  rotation=np.radians(label_position),ha='center',va='center')
ax1.grid(True)

# plot a quarter circle
circle1=plt.Circle((0,0),.99,color='k',transform=ax1.transData._b, alpha=0.7)
circle2=plt.Circle((0,0),.5,color='k',transform=ax1.transData._b)
ax1.add_artist(circle1)
ax1.add_artist(circle2)
#th = np.linspace(0,np.pi/2,100)
#ax1.plot(th, .99*np.ones_like(th), 'k-', linewidth=1)

#f1 = plt.figure()
#a1 = fractional_polar_axes(f1, thlim=(0,90), step=(15,0.2),
#  theta_offset=0)
#  
#a1.pcolormesh(theta,rad,data[0,:,:])
#ax = fig.add_subplot(111, projection="polar")

# TODO: clip polar plot
#p1 = ax.pcolormesh(X,Y,data[0,:,:])
#plt.colorbar(p1)
#ax.set_theta_direction(-1)
#ax.set_theta_zero_location("N")


imgname = imgdir + "part-pair-distri"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
