#!/usr/bin/env python2
import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt

# from floating_polar import fractional_polar_axes

print ""
print " ---- Particle Pair Distribution Plotting Utility ---- "
print ""

## Nicely sorted function
def sorted_nicely( l ):
  """ Sorts the given iterable in a natural way

  Required arguments:
  l -- The iterable to be sorted.

  courtesy stackoverflow/2669059
  """
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
  return sorted(l, key = alphanum_key)

# SIMULATION PARAMETERS
partR = 2.1

# Parse command line args and set up directory structure
if len(sys.argv) == 2:
  simdir = sys.argv[1]
else:
  print "Wrong number of inputs!"
  sys.exit()

if not simdir.endswith('/'):
  simdir = simdir + '/'

home = os.path.expanduser("~")
root = home + "/scratch/triply_per/" + simdir
datadir = root + "analysis/part-pair-distribution/3d-hist_bins/data/"

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

# Prompt for more input
nBinsTh = int(raw_input("Number of bins in theta: "))
nBinsR = int(raw_input("Number of bins in r: "))
compareFlag = int(raw_input("Greater Than (1) or Less Than (-1): "))
nSDEV = float(raw_input("Number of standard deviations: "))
combo = raw_input("Run a combo comparision? (y/n) ")

if compareFlag == 1:
  compare = "gtr"
elif compareFlag == -1:
  compare = "lss"
else:
  print "Compare flag needs to be (1) or (-1)"
  sys.exit()

# Data Files
infoFile = datadir + "bininfo-r%d-th%d" % (nBinsR, nBinsTh)

# Create imgdir if necessary
imgdir = root + "analysis/part-pair-distribution/3d-hist_bins/img"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Find time and evalZ, and size of each
time = np.genfromtxt(infoFile, skip_footer=2)[1:]
nt = np.size(time)

evalR = np.genfromtxt(infoFile, skip_header=1,skip_footer=1)[1:] /(2.*partR)
dr = np.mean(np.diff(evalR))
evalR = np.append(evalR, evalR[-1] + dr)

evalTh = np.genfromtxt(infoFile, skip_header=2)[1:]
dth = np.mean(np.diff(evalTh))
evalTh = np.append(evalTh, evalTh[-1] + dth)

## Either output one image, or compare a greater than / less than combo
if combo == "n":
  # -nBinsR-nBinsTh
  case = "bins-r%d-th%d_" % (nBinsR, nBinsTh)
  case = case + compare + "-%.1fsdev" % nSDEV
  data = np.genfromtxt(datadir + case)
  vmin = 0.
  vmax = 1.5
  print np.max(data)
  cmap = "viridis"
elif combo == "y":
  # compare a greater and less than case
  case = "bins-r%d-th%d_" % (nBinsR, nBinsTh)
  case_gtr = case + "gtr" + "-%.1fsdev" % np.abs(nSDEV)
  if (np.abs(nSDEV) < 1e-8):  # i.e. if 0
    case_less = case + "lss" + "-%.1fsdev" % np.abs(nSDEV)
  else:
    case_less = case + "lss" + "--%.1fsdev" % np.abs(nSDEV)
  data_gtr = np.genfromtxt(datadir + case_gtr)
  data_less = np.genfromtxt(datadir + case_less)
  data = data_gtr - data_less
  case = "compare-%.1fsdev" % nSDEV
  vmin = -np.max(np.abs(data));
  vmax = np.max(np.abs(data));
  cmap = "coolwarm"
else:
  print "Combo is" + combo
  sys.exit()

#data /= np.amax(data); # norm to 1

## Plot ###
theta,rad = np.meshgrid(evalTh, evalR)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, polar = True)
c = ax1.pcolormesh(theta, rad, data, norm=None,
  cmap=cmap, vmin=vmin, vmax=vmax)

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

# Set r-grids
rticks = np.array(ax1.get_yticks()).astype(int).astype(str)
rticks[::2] = ""
ax1.set_yticklabels(rticks)

# plot a quarter circle
circle1=plt.Circle((0,0),.99,color='k',transform=ax1.transData._b, alpha=0.7)
circle2=plt.Circle((0,0),.5,color='k',transform=ax1.transData._b)
ax1.add_artist(circle1)
ax1.add_artist(circle2)

imgname = imgdir + "/" + case 
print "      Printing to " + imgname
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

### Average over th ###
radial_distribution = np.mean(data, 1)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

evalR_centers = evalR[0:-1] + np.mean(np.diff(evalR))
ax2.plot(evalR_centers, radial_distribution, 'k')

ax2.set_xlim([0, evalR[-1]])
ax2.set_xlabel(r'$r/d$')
if combo == 0:
  ax2.set_ylim(ymin=0)
elif combo == 1:
  ax2.set_ylim([-np.max(np.abs(radial_distribution)), np.max(np.abs(radial_distribution))])
ax2.set_ylabel(r'$g(r)$')

imgname = imgdir + "/" + case + "_radial"
print "      Printing to " + imgname
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

