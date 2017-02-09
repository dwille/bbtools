#!/usr/bin/env python2

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
os.system('clear')


## Simulation Parameters
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

## Directories
home = os.path.expanduser("~")
root = home + "/scratch/collision-stokes_num/"
outfile = "bbottle.out"

## Init data arrays
st_mean = np.empty((rho.size, nparts.size))
st_mean[:,:] = np.nan
st_max = np.empty((rho.size, nparts.size))
st_max[:,:] = np.nan
st_sdev = np.empty((rho.size, nparts.size))
st_sdev[:,:] = np.nan
st_skew = np.empty((rho.size, nparts.size))
st_skew[:,:] = np.nan
st_num = np.empty((rho.size, nparts.size))
st_num[:,:] = np.nan

## Load data and analyze data -- loop over all output directories
for rr, curr_rho in enumerate(rho):
  for nn, curr_nparts in enumerate(nparts):
    simdir = root + str(curr_nparts) + "/rho" + str(curr_rho) + "/"

    # Proceed if simdir and outfile exist
    if os.path.exists(simdir) and os.path.exists(simdir + outfile):
      data = np.empty((0,3))

      # Open file and iterate through lines. If line indicates a collision, pull
      # out the part indices and Stokes numbers
      with open(simdir + outfile, "r") as f:
        for line in f.readlines():
          if line.startswith("COLL >> "):
            num = []
            for t in line.split():
              try:
                num.append(float(t))
              except ValueError:
                pass
            data = np.vstack((data, np.array(num)))

        # Sort first two columns of data and find unique contacts
        data[:,0:2] = np.sort(data[:,0:2], 1)
        data = np.vstack({tuple(row) for row in data})

        # Add data to arrays
        p_ind = np.argwhere(curr_rho == rho)
        n_ind = np.argwhere(curr_nparts == nparts)

        st_max[p_ind, n_ind] = data[:,2].max()
        st_mean[p_ind, n_ind] = data[:,2].mean()
        st_sdev[p_ind, n_ind] = data[:,2].std()
        st_skew[p_ind, n_ind] = stats.skew(data[:,2])
        st_num[p_ind, n_ind] = np.size(data,0)

print st_skew
print st_num

## Plotting
fig, ax = plt.subplots(1,1)

# mean vs. vfrac
for i in np.arange(rho.size):
  plt.plot(vfrac, st_mean[:,i], 'o--')

ax.legend([r"$\rho^* = 2.0$", r"$\rho^* = 3.3$", r"$\rho^* = 4.0$", 
  r"$\rho^* = 5.0$"], framealpha=0.7)

# mean + std vs vfrac
ax.set_prop_cycle(None)   # reset color cyace
for i in np.arange(rho.size):
  plt.plot(vfrac, st_mean[:,i] + st_sdev[:,i], ':')

# max
ax.set_prop_cycle(None)   # reset color cyace
for i in np.arange(rho.size):
  plt.plot(vfrac, st_max[:,i], '-.')

ax.set_xlabel(r"$\phi$")
ax.set_xlim([0, 0.5])
ax.set_ylabel(r"$St_c$")
ax.set_ylim(ymin=0)

# Save Figure
imgname = root + "/simdata/img/collision-stokes_num"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
