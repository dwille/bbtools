#!/usr/bin/env python2
import matplotlib.pyplot as plt
import numpy as np
import os
os.system('clear')

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
root = home + "/scratch/triply_per/"
datafile = "analysis/phase-averages/part/data/particleAvgVel"
stattime = "simdata/statTime"
stat_times = np.genfromtxt(root + stattime)

# Load data and find mean
for curr_nparts in nparts:
  for curr_rho in rho:
    simdir = root + str(curr_nparts) + "/rho" + str(curr_rho) + "/"

    # Find stat steady time
    is_n = np.argwhere(curr_nparts == stat_times[:,0])
    is_rho = np.argwhere(curr_rho == stat_times[:,1])
    i = np.intersect1d(is_n, is_rho)
    curr_stat_time = stat_times[i,2]

    # Proceed if simdir and outfile exist
    if os.path.exists(simdir) and os.path.exists(simdir + datafile):
      time = np.genfromtxt(simdir + datafile, skip_header=1, usecols=0)

      # Only take data after stat steady time
      t_start_ind = np.squeeze(np.argwhere(time >= curr_stat_time)[0])
      time = time[t_start_ind:]
      te = np.genfromtxt(simdir + datafile, skip_header=1, usecols=0)[t_start_ind:]
      wp = np.genfromtxt(simdir + datafile, skip_header=1, usecols=3)[t_start_ind:]

      print "%04d/rho%.1f:" % (curr_nparts, curr_rho)
      print "   time = [%.1f, %.1f] " % (time[0], time[-1])
      print "   <w_p> = %.5f" % wp.mean()



