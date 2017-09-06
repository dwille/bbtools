#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import matplotlib.lines as mlines
os.system('clear')

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

### Physical Parameters ###
d = 4.2                                                 # mm
nu = .01715                                             # mm^2/ms
g = 0.00981                                             # mm/ms^2
rho_f = 8.75e-4                                         # g/mm^3
rho_p = np.array([0.001750, 0.00289, 0.0035, 0.004375]) # g/mm^3
Lx = 42                                                 # mm
Ly = 42                                                 # mm
Lz = 126                                                # mm

### Hardcoded directories ###
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
imgdir = root + "img/f-rec-1D/"
freq_file = root + "wave-freq"

if not os.path.exists(freq_file):
  print("%s does not exist, exiting" % freq_file)
  sys.exit()

### Read and parse data ###
freq_data = np.genfromtxt(freq_file, skip_header=10, delimiter=",")

nparts = np.unique(freq_data[:,0]).astype(int)
rho = np.unique(freq_data[:,1])
phi = 4./3.*np.pi*(0.5*d)**3.*nparts/(Lx * Ly * Lz)

freq = np.zeros((rho.size, nparts.size))
for cc in np.arange(np.size(freq_data, 0)):
  curr_n = int(freq_data[cc, 0])
  curr_rho = freq_data[cc, 1]


  # Find indices to fill arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  freq[pp,nn] = freq_data[cc,2]

######################
### plotting #########
######################
size = (2,2)

### as function of phi ###
fig1 = plt.figure(figsize=size)
ax1 = fig1.add_subplot(111)

for i in np.arange(rho.size):
  plt.plot(phi, freq[i,:], 'o--')

ax1.set_xlabel(r"$\phi$")
ax1.set_xlim([0, 0.5])

ax1.set_ylabel(r"$f\ [Hz]$")
ax1.set_ylim([0, 2.25])

#ax1.legend([r"$\rho^* = 2.0$",r"$\rho^* = 3.3$",r"$\rho^* = 4.0$",r"$\rho^* = 5.0$"])

imgname = imgdir + "frequency_v_phi"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

### as function of rho ###
fig2 = plt.figure(figsize=size)
ax2 = fig2.add_subplot(111)

for i in np.arange(phi.size):
  plt.plot(rho, freq[:,i], 'o--')

ax2.set_xlabel(r"$\rho$")
ax2.set_xlim([1, 6])

ax2.set_ylabel(r"$f\ [Hz]$")
ax2.set_ylim([0, 2.25])

ax2.legend([r"$\rho^* = 2.0$",r"$\rho^* = 3.3$",r"$\rho^* = 4.0$",r"$\rho^* = 5.0$"])

imgname = imgdir + "frequency_v_rho"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
