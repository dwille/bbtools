#!/usr/bin/env python2
from setup import *

os.system('clear')

print ""
print " ---- Ternary Plot ---- "

# Setup directory structure and print
root = os.path.expanduser("~") + "/scratch/triply_per/"
imgdir = root + "simdata/img/tetrads/"
datadir = "/analysis/tetrads/data/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)
print " Root dir: " + root
print " Img dir:  " + imgdir

# Parameter sweep info
check_nparts = np.array([500, 1000, 1500, 2000])
check_rho = np.array([2., 3.3, 4., 5.])
simList = ['500/rho2.0', '500/rho3.3', '500/rho4.0', '500/rho5.0',
          '1000/rho2.0', '1000/rho3.3', '1000/rho4.0', '1000/rho5.0',
          '1500/rho2.0', '1500/rho3.3', '1500/rho4.0',
          '2000/rho2.0', '2000/rho3.3']
nSims = len(simList)

# Sort out colors
baseColors = ['r', 'g', 'b', 'k']
baseShades = [0.4, 0.57, 0.74, 0.9]
colors = ['']*nSims
shades = ['']*nSims
for cc, currsim in enumerate(simList):
  # Break directory string to grab nparts and rho
  nparts = int(currsim.split("/")[0])
  rho = float(currsim.split("/")[1][3:])

  # Different particle for different volume fractions
  for nn, n_check in enumerate(check_nparts):
    if (nparts == n_check):
      colors[cc] = baseColors[nn]

  # Different intensities for different density ratios
  for pp, p_check in enumerate(check_rho):
    if (rho == p_check):
      shades[cc] = baseShades[pp]

# Initialize data
nTetrads = np.zeros(nSims)
data = [ structtype() for i in range(nSims) ]
I1 = [ structtype() for i in range(nSims) ]
I2 = [ structtype() for i in range(nSims) ]
I3 = [ structtype() for i in range(nSims) ]

# Loop over all directory's
legendText = ['']*nSims
for ii in np.arange(nSims):
  simdir = root + simList[ii] + datadir
  infoFile = simdir + 'info.dat'
  nodeFile = simdir + 'regularNodes'
  statmean = simdir + 'stat.mean'

  nTetrads[ii] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

  # Pull time [ms]
  tmpT = np.genfromtxt(statmean, skip_header=1, usecols=0)
  data[ii].time = np.zeros(np.size(tmpT))
  data[ii].time = tmpT - tmpT[0]

  # Init arrays
  I1[ii].mean    = np.zeros(np.size(tmpT))
  I2[ii].mean    = np.zeros(np.size(tmpT))
  I3[ii].mean    = np.zeros(np.size(tmpT))

  # Pull data
  I1[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=4)
  I2[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=5)
  I3[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=6)

  legendText[ii] = simList[ii] # + ': ' + str(nTetrads[ii])

# Ternary plot
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)

# Geometry
L = 1.
H = np.sqrt(3.)/2.

# Limits and ticks
ax.set_xlim([0,L])
ax.set_xticks([])
ax.set_ylim([0,L])
ax.set_yticks([])
ax.axis('off')

# Outer lines
ax.plot([0, L], [0, 0], 'k-', linewidth=4)
ax.plot([0, 0.5*L], [0, H], 'k-', linewidth=2)
ax.plot([0.5*L, L], [H, 0], 'k-', linewidth=2)

# Axes
ax.plot([0, H*H], [0, 0.5*H], 'k--')
ax.plot([L, 1. - H*H], [0., 0.5*H], 'k--')
ax.plot([0.5*L, 0.5*L], [0, H], 'k--')

# Center
ax.plot(0.5*L, H/3., 'ko')

# Labels
plt.text(0, -0.06,             r"$I_1 = 1$", color='blue')
plt.text(0.78, np.sqrt(3.)/4., r"$I_1 = 0$", color='blue')
plt.text(1., -0.06,            r'$I_2 = 1$', color='red')
plt.text(0.07, np.sqrt(3.)/4., r'$I_2 = 0$', color='red')
plt.text(0.5, 0.9,             r'$I_3 = 1$', color='green')
plt.text(0.5, -.06,            r'$I_3 = 0$', color='green')

# Grid lines? (see bbtools/deprecated/matlab/devel/ternary_fig.m)

# Actually plot...
for ii in np.arange(nSims):
  # Convert to correct coords
  a = I1[ii].mean
  b = I2[ii].mean
  c = I3[ii].mean
  xx = 0.5*(2.*b + c)/(a + b + c);
  yy = 0.5*np.sqrt(3.)*c/(a + b + c);

  ax.plot(xx, yy, color=colors[ii], alpha=shades[ii], linewidth=.5)


# Save
imgname = imgdir + "ternary_inorm"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
