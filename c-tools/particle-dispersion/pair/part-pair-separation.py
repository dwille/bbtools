#!/usr/bin/env python

# PURPOSE
#   Calculate the autocorrelation of particle velocity from the supplied time
#   series.
# USAGE
#   ./part-pair-separation.py </path/to/sim/output> <start_time>
#
# OUTPUT
#   </path/to/sim/analysis/...XXX

# Imports:
import sys
import numpy as np
from scipy import signal
import scipy.spatial as spatial
import bluebottle_particle_reader as bbparts

##########

# Parse output directory from command line
if len(sys.argv) >= 2:    # output directory given
  data_dir = sys.argv[1]

  if len(sys.argv) >= 3:  # start time given
    t_start = float(sys.argv[2])

else:                     # nothing given
  print("Error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./part-pair-separation.py <./path/to/sim/output> <start_time>")
  print(" or")
  print("   ./part-pair-separation.py <./path/to/sim/output>")
  sys.exit()

# Initialize the reader, grab times in range
times = bbparts.init(data_dir)
time_flt = np.array([tt for tt in times], dtype='float')
ind = np.argwhere(time_flt >= t_start)[0][0]
times = times[ind:]
print("  Running from t = %s to t = %s" % (times[0], times[-1]))
time_flt = time_flt[ind:] - time_flt[ind]

# Get number of particles
bbparts.open(times[0])
nparts = bbparts.read_nparts()
bbparts.close()

# Predefine simulation parameters -- XXX should not hardcode!
a = 2.1
Lx = 42.
Ly = 42.
Lz = 126.
r = 3.*a

###### Initialize particle pairs using a k-d-tree ######

# Pull initial particle positions
bbparts.open(times[0])
(x0, y0, z0) = bbparts.read_part_position()
bbparts.close()

# If we initialize over periodic boundaries, then flipping later for periodic
# boundaries is difficult. We did it for tetrads, but for now, just ignore that.
# Rescale data to ([0,Lx], [0, Ly], [0, Lz])
#boxsize = np.array([Lx,Ly,Lz,Lx,Ly,Lz])
#data = np.c_[x0.ravel() + 0.5*Lx,
#             y0.ravel() + 0.5*Ly,
#             z0.ravel() + 0.5*Lz]
data = np.c_[x0.ravel(), y0.ravel(), z0.ravel()]

# Set up tree with particle positions
tree = spatial.cKDTree(data)#, boxsize=boxsize)

# Query particle pairs less than desired radius
pairs = tree.query_pairs(r)
npairs = len(pairs)
print("  Found %d pairs within a distance of %.2f" % (npairs, r))

# Find initial separation -- could also use cKDTree.spase_distance_matrix here
init_distance = np.zeros(npairs)
for pp,pair in enumerate(pairs):
  i = pair[0]
  j = pair[1]

  dx = x0[i] - x0[j]
  #dx += Lx * (1*(dx < -r) - 1*(dx > r))  # if init over boundaries

  dy = y0[i] - y0[j]
  #dy += Ly * (1*(dy < -r) - 1*(dy > r))  # if init over boundaries

  dz = z0[i] - z0[j]
  #dz += Lz * (1*(dz < -r) - 1*(dz > r))  # if init over boundaries

  init_distance[pp] = np.sqrt(dx**2 + dy**2 + dz**2)

print("  Mean initial separation: %.2f" % np.mean(init_distance))
print("  Sdev initial separation: %.2f" % np.std(init_distance))
print("  Min initial separation:  %.2f" % np.min(init_distance))

###### Loop over time and find mean squared separation ######
r2_total = np.zeros(len(times))
r2_horiz = np.zeros(len(times))
r2_verti = np.zeros(len(times))
theta_z_mean = np.zeros(len(times))
theta_z_sdev = np.zeros(len(times))
theta_z_max = np.zeros(len(times))
theta_z_min = np.zeros(len(times))

flip_i = np.zeros(nparts)
flip_j = np.zeros(nparts)
flip_k = np.zeros(nparts)

dz_time_trace = np.zeros((len(times), npairs))

for tt,time in enumerate(times):
  print("  Time %d of %d\r" % (tt+1, len(times)), end="")
  bbparts.open(time)
  (x, y, z) = bbparts.read_part_position();

  ### Deal with particles crossing periodic boundaries ###
  # See note above about initializing across boundaries. For now, don't do that
  # from single particle fun_kernel

  # Before anything, correct particle position with old flip count to make sure
  # we're comparing apples to apples
  x += Lx * flip_i
  y += Ly * flip_j
  z += Lz * flip_k

  # Position change since previous timestep
  dx = x0 - x
  dy = y0 - y
  dz = z0 - z

  # Increment flip count:
  #  -- if (prev - curr) > 0, went R->L, (+) a domain length
  #  -- if (prev - curr) < 0, went L->R, (-) a domain length
  flip_i += 1*(dx >= 0.5*Lx) - 1*(dx <= -0.5*Lx)
  flip_j += 1*(dy >= 0.5*Ly) - 1*(dy <= -0.5*Ly)
  flip_k += 1*(dz >= 0.5*Lz) - 1*(dz <= -0.5*Lz)

  # Correct particle position one last time
  x += Lx * (1*(dx >= 0.5*Lx) - 1*(dx <= -0.5*Lx))
  y += Ly * (1*(dy >= 0.5*Ly) - 1*(dy <= -0.5*Ly))
  z += Lz * (1*(dz >= 0.5*Lz) - 1*(dz <= -0.5*Lz))

  # Loop over pairs and find separation
  # XXX better way -- vectorize?
  dx_pair = np.zeros(npairs)
  dy_pair = np.zeros(npairs) 
  dz_pair = np.zeros(npairs) 
  for pp,pair in enumerate(pairs):
    i = pair[0]
    j = pair[1]

    dx_pair[pp] = x[i] - x[j]
    dy_pair[pp] = y[i] - y[j]
    dz_pair[pp] = z[i] - z[j]

  r2_total[tt] = np.mean(dx_pair**2 + dy_pair**2 + dz_pair**2)
  r2_verti[tt] = np.mean(dz_pair**2)
  r2_horiz[tt] = np.mean(dx_pair**2 + dy_pair**2)

  # Find alignment with vertical
  # cos(theta_z) = (r_ab,e_z) / ||r_ab|| ||e_z|| = dz / ||r_ab||
  # Not concerned with up/down, so do |dz|
  arg = np.abs(dz_pair) / np.sqrt(dx_pair**2 + dy_pair**2 + dz_pair**2)
  theta_z_mean[tt] = np.mean(arg)
  theta_z_sdev[tt] = np.std(arg)
  theta_z_max[tt] = np.max(arg)
  theta_z_min[tt] = np.min(arg)

  #
  dz_time_trace[tt,:] = dz_pair
 
  bbparts.close()

  # Save data as old timestep
  x0 = x
  y0 = y
  z0 = z

# XXX alignment
# XXX add plus or minus sdev to plot?
 
print()

# Temporary plot
import matplotlib.pyplot as plt
plt.figure()

a = 2.1 # mm
rho = 3.3
nu = 0.01715
tau_p = (2.*a)**2 * rho / (18. * nu)

plt.loglog(time_flt/tau_p, np.sqrt(r2_total)/a, '.', markersize=1)
plt.loglog(time_flt/tau_p, np.sqrt(r2_verti)/a, '.', markersize=1)
plt.loglog(time_flt/tau_p, np.sqrt(r2_horiz)/a, '.', markersize=1)
plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$\sqrt{\langle r_{\alpha \beta}^2\rangle}/a$")
plt.legend(["total", "verti", "horiz"])
plt.tight_layout()
plt.savefig("tmp-fig.png", bbox_inches='tight', format='png')

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.semilogx(time_flt/tau_p, theta_z_mean, '.', markersize=1)
ax1.semilogx(time_flt/tau_p, theta_z_mean + theta_z_sdev, 'k.', markersize=1)
ax1.semilogx(time_flt/tau_p, theta_z_mean - theta_z_sdev, 'k.', markersize=1)
ax1.semilogx(time_flt/tau_p, theta_z_max, 'r.', markersize=1)
ax1.semilogx(time_flt/tau_p, theta_z_min, 'r.', markersize=1)
ax1.set_xlabel(r"$t/\tau_p$")
ax1.set_ylabel(r"$\langle \cos(\theta_z) \rangle_{n_p}$")
ax1.set_ylim([0, 1])

# To plot axis on right side with theta not cos(theta)
#ax2 = ax1.twinx()
#ax2.semilogx(time_flt/tau_p, np.arccos(theta_z), '.', markersize=0)
#ax2.set_ylabel(r"$\theta_z$")
#ax2.set_ylim([0, 0.5*np.pi])

plt.tight_layout()
plt.savefig("tmp-fig-ang.png", bbox_inches='tight', format='png')

# 
fig = plt.figure()
for pp,_ in enumerate(pairs):
  plt.plot(times, dz_time_trace[:,pp])

plt.xlabel("time")
plt.ylabel("dz")
plt.savefig("tmp-dz-trace.png")
