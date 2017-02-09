#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import matplotlib.lines as mlines
os.system('clear')

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

print ""
print " ---- Wavespeed Plotting Utility ---- "
print ""

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

# Read data
wave_data = np.genfromtxt(root + "waveData_normalized", skip_header=2)

nparts = np.unique(wave_data[:,0])
rho = np.unique(wave_data[:,1])
phi = 4./3.*np.pi*(0.5*d)**3.*nparts/(Lx * Ly * Lz)

# Pull wavespeeds to convenient arrays
c_phi = np.zeros((rho.size, nparts.size))
c_wp = np.zeros((rho.size, nparts.size))
for cc in np.arange(np.size(wave_data, 0)):
  curr_n = wave_data[cc, 0]
  curr_rho = wave_data[cc, 1]

  # Find indices to fill arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  c_phi[pp,nn] = wave_data[cc,2]
  c_wp[pp,nn] = wave_data[cc,3]

## Wave Speed Plot ##
fig1 = plt.figure(figsize=(3.25,3.25))
ax1 = fig1.add_subplot(111)

for i in np.arange(rho.size):
  plt.plot(phi, c_wp[i,:]/c_phi[i,:], 'o--')

ax1.set_xlabel(r'$\phi$')
ax1.set_xlim([0, 0.5])

ax1.set_ylabel(r'$c_{w_p}^* / c_{\phi}^*$')
ax1.set_ylim([1, 1.25])

ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_major_locator(MultipleLocator(.05))
ax1.yaxis.set_minor_locator(MultipleLocator(.025))

ax1.legend([r"$\rho^* = 2.0$",r"$\rho^* = 3.3$",r"$\rho^* = 4.0$",r"$\rho^* = 5.0$"])

imgname = imgdir + "wavespeed_comp_phi_wp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
#plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')
