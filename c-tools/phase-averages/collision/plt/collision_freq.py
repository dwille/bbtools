#!/usr/bin/env python2
import sys, os, glob, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

## Struct class
class structtype():
  pass

# Directory structure
home = os.path.expanduser("~")
root = home + "/scratch/collision/"
datadir = root + "/*/rho*/analysis/phase-averages/collision/data/collision_stats"
imgdir = root + "/simdata/img/colls/"

if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Sims that we have data for
simdirs = glob.glob(datadir)
nsims = int(len(simdirs))

# Pull data
data = [ structtype() for i in range(nsims) ]
for cc, currdir in enumerate(simdirs):
  # Pull data, cut out first time step
  data[cc].time = np.genfromtxt(currdir, skip_header=1, usecols=0)
  data[cc].ncolls = np.genfromtxt(currdir, skip_header=1, usecols=1)

  nsteps = np.size(data[cc].time)

  freq = np.zeros(nsteps)
  dt = np.zeros(nsteps)
  for ii in np.arange(1,nsteps):
    dcoll = data[cc].ncolls[ii] - data[cc].ncolls[0]
    dt[ii] = data[cc].time[ii] - data[cc].time[0]
    freq[ii] = dcoll/dt[ii]

  data[cc].freq = freq
  data[cc].dt = dt

# colors for plotting -- still ungeneral but better than in tetrads
check_nparts = np.array([500,1000,1500,2000])
check_rho = np.array([2., 3.3, 4., 5.])
baseColors = ['r', 'g', 'b', 'k']
baseShades = [0.4, 0.57, 0.74, 0.9]
colors = ['']*nsims
shades = ['']*nsims
for cc, currdir in enumerate(simdirs):
  # Break directory string to grab nparts and rho
  nparts = int(currdir.split("/")[5])
  rho = float(currdir.split("/")[6][3:])

  # Different particle for different volume fractions
  for nn, n_check in enumerate(check_nparts):
    if (nparts == n_check):
      colors[cc] = baseColors[nn]

  # Different intensities for different density ratios
  for pp, p_check in enumerate(check_rho):
    if (rho == p_check):
      shades[cc] = baseShades[pp]

# plot
plt.figure()

for cc in np.arange(nsims):
  plt.plot(data[cc].dt, data[cc].freq, color=colors[cc], alpha=shades[cc])
  #plt.plot(data[cc].time - data[cc].time[0], 0.5*(data[cc].ncolls-data[cc].ncolls[0]),
  # color=colors[cc], alpha=shades[cc], marker='o')

plt.xlabel(r"$t$")
plt.ylabel(r"$coll/t$")

# save
imgname = imgdir + "coll_freq"
print "Saving figure to %s" % imgname
plt.savefig(imgname + ".png", bbox_inches="tight", format='png')
