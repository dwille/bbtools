#!/usr/bin/env python2

import matplotlib
matplotlib.use(u'Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import os,sys
#from scipy.misc import imread
from matplotlib.image import imread

# Simulation Parameters
d = 4.2                                                 # mm
nu = .01715                                             # mm^2/ms
g = 0.00981                                             # mm/ms^2
rho_p = np.array([0.001750, .00289, 0.0035, 0.004375])  # g/mm^3
rho_f = 8.75e-4                                         # mm
vfrac = np.array([0.087, 0.175, 0.262, 0.349])          # []
rho = np.array([2.0, 3.3, 4.0, 5.0])                    # []
nparts = np.array([500, 1000, 1500, 2000])              # []
Lx = 42                                                 # mm
Ly = 42                                                 # mm
Lz = 126                                                # mm

# Directories
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
imgdir = root + "img/velocity_flucts/"

# Load fluctuations data
# u',v',w'       -- standard deviation of particle velocities [mm/ms]
# urel,vrel,wrel -- calculated from phase averaged velocities [mm/ms]
fluct_data = np.genfromtxt(root + "velFlucts_raw", skip_header=1)

# Calculate Stokes' velocity -- mm/ms
stokes_vel = 2./9.*(rho - 1.)*(0.5*d)*(0.5*d)*g/nu
stokes_Reynolds = stokes_vel * 0.5*d / nu

# Pull data to arrays
u = np.zeros((rho.size, nparts.size))
v = np.zeros((rho.size, nparts.size))
w = np.zeros((rho.size, nparts.size))
for cc in np.arange(np.size(fluct_data,0)):
  curr_n = fluct_data[cc, 0]
  curr_rho = fluct_data[cc, 1]

  # Find indices to fill arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  # Normalize fluctuations by relative velocity
  mean_relative_velocity = np.sqrt(fluct_data[cc,5]**2. 
                                 + fluct_data[cc,6]**2. 
                                 + fluct_data[cc,7]**2.)
  u[pp,nn] = fluct_data[cc, 2] / mean_relative_velocity
  v[pp,nn] = fluct_data[cc, 3] / mean_relative_velocity
  w[pp,nn] = fluct_data[cc, 4] / mean_relative_velocity

  #print "%d -- %.1lf: %.3lf" % (curr_n, curr_rho, stokes_vel[pp]*stokes_Reynolds[pp]**(-1./3.)/mean_relative_velocity)
  #print "%d -- %.1lf: %.3lf" % (curr_n, curr_rho, 2.*stokes_vel[pp]/mean_relative_velocity)

# Plots #
fig1 = plt.figure(figsize=(5,2))

# velocity vs vfrac
ax1 = fig1.add_subplot(121)
for i in np.arange(rho.size):
  plt.loglog(vfrac, w[i,:], 'o--')

ax1.set_xlabel(r"$\phi$")
ax1.set_xlim([0.05, 0.5])
ax1.set_ylabel(r"$w/\langle w \rangle$")
#ax1.set_ylabel(r"$w$")
ax1.set_ylim(0.1, 1)
ax1.set_title(r"$w$")
ax1.xaxis.set_major_locator(MultipleLocator(0.1))

xpts = [0.06, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
ypts = 0.5*np.power(xpts, 1./3.)
plt.loglog(xpts, ypts, 'k--')

ax2 = fig1.add_subplot(122)
ax2.set_prop_cycle(None)
for i in np.arange(rho.size):
  #plt.semilogx(vfrac, u[i,:], 'o--')
  plt.loglog(vfrac, u[i,:], 'o--')

ax2.set_prop_cycle(None)
for i in np.arange(rho.size):
  #plt.semilogx(vfrac, v[i,:], 's:')
  plt.plot(vfrac, v[i,:], 's:')

ax2.set_xlabel(r"$\phi$")
ax2.set_xlim([0.05, 0.5])
ax2.set_title(r"$u,v$")
ax2.yaxis.set_ticklabels([])
ax2.set_ylim(0.1, 1)
ax2.xaxis.set_major_locator(MultipleLocator(0.1))
#ax2.set_ylabel(r"${u,v}^\prime/\langle w \rangle$")

# Save fig
#imgname = imgdir + "velocity_flucts_linear"
#imgname = imgdir + "velocity_flucts_linear_unnormalized"
imgname = imgdir + "velocity_flucts_log"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

##############################
# Find and load histogram data
hist_data = np.genfromtxt(root + "vfrac_vel_hist", skip_header=2)

# Find number of particles, density ratios -> volume fraction
nparts = np.unique(hist_data[:,0])
rho = np.unique(hist_data[:,1])
vfrac = 4./3.*np.pi*(0.5*d)**3.*nparts/(Lx * Ly * Lz)

# Pull slope and r values to convenient arrays
m = np.zeros((rho.size, nparts.size))
for cc in np.arange(np.size(hist_data, 0)):
  curr_n = hist_data[cc, 0]
  curr_rho = hist_data[cc, 1]

  # Find indices to fill m,r arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  mean_relative_velocity = np.sqrt(fluct_data[cc,5]**2. 
                                 + fluct_data[cc,6]**2. 
                                 + fluct_data[cc,7]**2.)
  m[pp,nn] = hist_data[cc, 5]

## Plots ##
fig1 = plt.figure(figsize=(2,2))

# slope vs vfrac
ax1 = fig1.add_subplot(111)
for i in np.arange(rho.size):
  plt.loglog(vfrac, m[i,:], 'o--')

xpts = [.07, .4]
ypts = 0.14*np.power(xpts, -2./3.)
plt.loglog(xpts, ypts, 'k--')

ax1.set_xlabel(r'$\phi$')
ax1.set_xlim([.05, 0.5])
ax1.set_ylabel(r'$dw/d\phi$')
ax1.set_ylim([.1,1])
#ax1.yaxis.set_major_locator(MultipleLocator(0.25))

ax1.arrow(0.1, 0.25, 0.06, 0.55, head_width=0.02, head_length=0.1, fc='k', ec='k')
ax1.text(0.2, 0.8, r"$\rho^*\uparrow$")

# Save fig
imgname = imgdir + "vfrac_vel_hist_flucts"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# ## Data overlay ##
# overlay_img = imgdir + "2011-hinch-guaz_vel_flucts.png"
# #overlay_img = imgdir + "2011-guazzelli_hinch-vel_flucts-cropped.png"
# img = imread(overlay_img)
# fig3 = plt.figure(figsize=(2,2))
# 
# for i in np.arange(rho.size):
#   plt.semilogx(vfrac, w[i,:], 'o--')
# plt.xlim([0.0001, 0.6])
# 
# plt.imshow(img, zorder=1, aspect="auto", interpolation="none", extent=[0.0001, 0.6, 0, 2])
# 
# # Save fig
# imgname = imgdir + "vel_flucts_overlay"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
# #plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Plot velocity fluctuations against wavespeed
wave_data = np.genfromtxt(root + "waveData", skip_header=1) ## mm/s

# Pull data to arrays
c = np.empty((rho.size, nparts.size))
c.fill(np.nan)
for cc in np.arange(np.size(wave_data,0)):
  curr_n = wave_data[cc, 0]
  curr_rho = wave_data[cc, 1]

  # Find indices to fill arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  # Normalize fluctuations by relative velocity
  tmp_c = wave_data[cc, 2] 
  if (tmp_c != -1):
    tmp_c /= 1000.       # mm/s -> mm/ms
    tmp_c *= d/nu # normalize by nu/d
    c[pp,nn] = tmp_c

fig2 = plt.figure(figsize=(2,2))
ax1 = fig2.add_subplot(111)

for i in np.arange(rho.size):
  ax1.plot(w[i,:], c[i,:], 'o--')

ax1.arrow(0.30, 11, -.05, 20, head_width=0.01, head_length=1, fc='k', 
  ec='k', linewidth=0.5)
ax1.text(0.29, 15, r"$\rho^*\uparrow$")

ax1.arrow(0.25, 5, .15, 5, head_width=0.03, head_length=.4, fc='k', 
  ec='k', linewidth=0.5)
ax1.text(0.35, 5, r"$\phi\uparrow$")

ax1.set_xlabel(r"$w'/\langle w \rangle$")
ax1.set_ylabel(r"$2ac/\nu$")
ax1.set_ylim([0,45])


# Save fig
imgname = imgdir + "velflucts_wavespeed"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')



