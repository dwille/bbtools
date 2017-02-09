#!/usr/bin/env python2
from setup import *
from matplotlib.ticker import MultipleLocator
os.system('clear')

# Input parameters
d = 4.2                                                 # mm
nu = .01715                                             # mm^2/ms
g = 0.00981                                             # mm/ms^2
rho_f = 8.75e-4                                         # g/mm^3
rho_p = np.array([0.001750, 0.00289, 0.0035, 0.004375]) # g/mm^3
Lx = 42                                                 # mm
Ly = 42                                                 # mm
Lz = 126                                                # mm

# Directory structure
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
imgdir = root + "img/f-rec-1D/"

# Find and load histogram data
hist_data = np.genfromtxt(root + "vfrac_vel_hist", skip_header=2)

# Find number of particles, density ratios -> volume fraction
nparts = np.unique(hist_data[:,0])
rho = np.unique(hist_data[:,1])
vfrac = 4./3.*np.pi*(0.5*d)**3.*nparts/(Lx * Ly * Lz)

# Pull slope and r values to convenient arrays
m = np.zeros((rho.size, nparts.size))
r = np.zeros((rho.size, nparts.size))
for cc in np.arange(np.size(hist_data, 0)):
  curr_n = hist_data[cc, 0]
  curr_rho = hist_data[cc, 1]

  # Find indices to fill m,r arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  m[pp,nn] = hist_data[cc, 5]
  r[pp,nn] = hist_data[cc,4]

## Plots ##
fig1 = plt.figure(figsize=(4,1.75))

# slope vs vfrac
ax1 = fig1.add_subplot(121)
for i in np.arange(rho.size):
  plt.loglog(vfrac, m[i,:], 'o--')

xpts = [.07, .4]
ypts = 0.14*np.power(xpts, -2./3.)
plt.loglog(xpts, ypts, 'k--')

ax1.set_xlabel(r'$\phi$')
ax1.set_xlim([.05, 0.5])
#ax1.xaxis.set_ticklabels([])
ax1.set_ylabel(r'$dw/d\phi$')
ax1.set_ylim([.1,1.])
#ax1.yaxis.set_major_locator(MultipleLocator(0.25))

ax1.arrow(0.1, 0.25, 0.06, 0.55, head_width=0.02, head_length=0.1, fc='k', ec='k')
ax1.text(0.2, 0.8, r"$\rho^*\uparrow$")

# slope vs density ratio
ax2 = fig1.add_subplot(122)
for i in np.arange(vfrac.size):
  plt.loglog(rho, m[:,i], 'o--')

ax2.set_xlabel(r"$\rho^*$")
ax2.set_xlim([1, 10])
#ax2.xaxis.set_ticklabels([])
ax2.yaxis.set_ticklabels([])
ax2.set_ylim([.1,1])
#ax2.yaxis.set_major_locator(MultipleLocator(0.25))

ax2.arrow(2.5, 0.70, 0.0, -0.55, head_width=0.15, head_length=0.04, fc='k', ec='k')
ax2.text(2.8, 0.13, r"$\phi\uparrow$")

# # r vs vfrac
# ax3 = fig1.add_subplot(223)
# for i in np.arange(rho.size):
#   plt.plot(vfrac, r[i,:], 'o--')
# 
# ax3.set_xlabel(r'$\phi$')
# ax3.set_xlim([0, 0.5])
# ax3.set_ylabel(r'$r$')
# ax3.set_ylim([0.75,1])
# ax3.yaxis.set_major_locator(MultipleLocator(0.05))
# 
# ax3.arrow(0.13, 0.76, -.03, 0.1, head_width=0.02, head_length=0.04, fc='k', ec='k')
# ax3.text(0.12, 0.9, r"$\rho^*\uparrow$")
# 
# # r vs density ratio
# ax4 = fig1.add_subplot(224)
# for i in np.arange(vfrac.size):
#   plt.plot(rho, r[:,i], 'o--')
# 
# ax4.set_xlabel(r"$\rho^*$")
# ax4.set_xlim([1, 6])
# ax4.yaxis.set_ticklabels([])
# ax4.set_ylim([0.75,1])
# ax4.yaxis.set_major_locator(MultipleLocator(0.05))
# 
# ax4.arrow(2.5, 0.76, 0, 0.16, head_width=0.15, head_length=0.02, fc='k', ec='k')
# ax4.text(2.8, 0.95, r"$\phi\uparrow$")

# Save fig
imgname = imgdir + "vfrac_vel_hist_results"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

