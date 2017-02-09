#!/usr/bin/env python2
from setup import *

# Sims that we have data for
simdirs = glob.glob(datadir)
nsims = int(len(simdirs))

# Data arrays -- [rho, nparts]
data = np.empty([rho_star.size, nparts.size])
data.fill(np.nan)

# Pull data
for cc, currdir in enumerate(simdirs):
  # Pull data, cut out first time step
  time = np.genfromtxt(currdir, skip_header=1, usecols=0)[1:]
  ncolls = np.genfromtxt(currdir, skip_header=1, usecols=1)[1:]
  freq =(ncolls[-1] - ncolls[0]) / (time[-1] - time[0]) 

  curr_nparts = int(currdir.split("/")[5])
  curr_rho = float(currdir.split("/")[6][3:])

  # Find an index to store in data array and store
  pp = np.argwhere(curr_rho == rho_star)
  nn = np.argwhere(curr_nparts == nparts)

  data[pp, nn] = freq

# colors for plotting -- still ungeneral but better than in tetrads
baseColors = ['r', 'g', 'b', 'k']
baseShades = [0.4, 0.57, 0.74, 0.9]
colors = ['']*nsims
shades = ['']*nsims
for cc, currdir in enumerate(simdirs):
  # Break directory string to grab nparts and rho
  curr_nparts = int(currdir.split("/")[5])
  curr_rho = float(currdir.split("/")[6][3:])

  # Different color for different volume fractions
  for nn, n_check in enumerate(nparts):
    if (curr_nparts == n_check):
      colors[cc] = baseColors[nn]

  # Different intensities for different density ratios
  for pp, p_check in enumerate(rho_star):
    if (curr_rho == p_check):
      shades[cc] = baseShades[pp]

# plot
fig1 = plt.figure(figsize=(4,6))

# Constant volume fraction, changing density ratio
ax1 = fig1.add_subplot(211)
plt.plot(rho_star, data[:,0], 'o--')
plt.plot(rho_star, data[:,1], 'o--')
plt.plot(rho_star, data[:,2], 'o--')
plt.plot(rho_star, data[:,3], 'o--')

plt.legend([r"$\phi = 0.087$",r"$\phi = 0.175$",r"$\phi = 0.262$",r"$\phi = 0.349$"],
  loc="upper right", framealpha=0.6)
plt.xlabel(r"$\rho^*$")
plt.ylabel(r"collisional frequency, $n_{coll}/ms$")
plt.xlim([1,6])

## Constant density ratio, changing volume fraction
ax2 = fig1.add_subplot(212)
plt.loglog(vfrac, data[0,:], 'o--')
plt.loglog(vfrac, data[1,:], 'o--')
plt.loglog(vfrac, data[2,:], 'o--')
plt.loglog(vfrac, data[3,:], 'o--')

plt.legend([r"$\rho^* = 2.0$",r"$\rho^* = 3.3$",r"$\rho^* = 4.0$",r"$\rho^* = 5.0$"],
  loc="lower right")
plt.xlabel(r"$\phi$")
plt.ylabel(r"collisional frequency, $n_{coll}/ms$")
plt.xlim([.05,1])
plt.ylim(ymax=125)

xpts = [.07, .50]
ypts = 1500*np.power(xpts, 3.)
print xpts
print ypts
plt.plot(xpts, ypts, 'k--')
plt.text(.07, .3, r"slope=3")

# save
imgname = imgdir + "coll_freq"
print "Saving figure to %s" % imgname
plt.savefig(imgname + ".png", bbox_inches="tight", format='png')
