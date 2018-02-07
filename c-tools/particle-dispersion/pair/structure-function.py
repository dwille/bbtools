#!/usr/bin/env python

# PURPOSE
#   Calculate the autocorrelation of particle velocity from the supplied time
#   series.
# USAGE
#   ./structure-function.py </path/to/sim/output> <start_time>
#
# OUTPUT
#   </path/to/sim/analysis/...XXX

# Imports:
import sys
import numpy as np
import scipy.spatial as spatial
import numpy.linalg as linalg

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
  print("   ./structure-function.py <./path/to/sim/output> <start_time>")
  print(" or")
  print("   ./structure-function.py <./path/to/sim/output>")
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
r = 20.*a

# Pull initial particle positions
bbparts.open(times[0])
(x0, y0, z0) = bbparts.read_part_position()
bbparts.close()

####################################
#### Use kd-tree to find pairs  ####
####################################
#XXX should we re-init particle each timestep or track them?
# If we initialize over periodic boundaries, then flipping later for periodic
# boundaries is difficult. We did it for tetrads, but for now, just ignore that.
# Rescale data to ([0,Lx], [0, Ly], [0, Lz])
#boxsize = np.array([Lx,Ly,Lz,Lx,Ly,Lz])
#data = np.c_[x0.ravel() + 0.5*Lx,
#             y0.ravel() + 0.5*Ly,
#             z0.ravel() + 0.5*Lz]

# Set up tree with particle positions
data = np.c_[x0.ravel(), y0.ravel(), z0.ravel()]
tree = spatial.cKDTree(data)#, boxsize=boxsize)

# Query particle pairs less than desired radius
pairs = tree.query_pairs(r)
pairs = np.array(list(pairs))
npairs = len(pairs)
print("  Found %d pairs within a distance of %.2f" % (npairs, r))

####################################################
#### Find and bin part separations and velocity ####
####################################################
### XXX NOTE THESE ARE EULERIAN, WE SHOULD DO LAGRANGIAN

# Bin parameters
dr = 0.25*a           # Distance between bin centers
bin_width = 0.5*a    # Width of bin. Should be >= dr,
bin_center = np.arange(2.*a + 0.5*dr, 10.*r, dr)
nbins = len(bin_center)

# These will hold the sum of increments -- divide by counts later to find mean
du_verti = np.zeros(nbins)            # Velocity increments in vertical
du_horiz = np.zeros(nbins)            # Velocity increments in horizontal
du_parallel = np.zeros(nbins)         # Vel increments parallel to sep
du_perpendicular = np.zeros(nbins)    # Vel increments perpendicular to sep

# These hold bin counts
du_bin_counts = np.zeros(nbins)

# Periodic flip counters
flip_i = np.zeros(nparts)
flip_j = np.zeros(nparts)
flip_k = np.zeros(nparts)

# Loop and find position and velocity at each timestep
for tt,time in enumerate(times):
  print("  Time %d of %d\r" % (tt+1, len(times)), end="")
  bbparts.open(time)
  (x, y, z) = bbparts.read_part_position();
  (u, v, w) = bbparts.read_part_velocity();

  # Deal with particles crossing periodic boundaries
  (x,y,z,flip_i,flip_j,flip_k) = bbparts.periodic_flip(x, y, z, x0, y0, z0,    \
                                                        flip_i, flip_j, flip_k,\
                                                        Lx, Ly, Lz)

  # Find position and separation of each pair
  r1 = np.array([x[pairs[:,0]], y[pairs[:,0]], z[pairs[:,0]]]).T
  r2 = np.array([x[pairs[:,1]], y[pairs[:,1]], z[pairs[:,1]]]).T
  r12 = r2 - r1

  # Sanity check using 3^2 + 4^2 + 12^2 = 13^2
  #r1[:,0] = 0; r1[:,1] = 0; r1[:,2] = 0;
  #r2[:,0] = 3; r2[:,1] = 4; r2[:,2] = 12;
  #r12 = r2 - r1

  r_mag = linalg.norm(r12, axis=1)

  # Find velocity difference for each pair
  u1 = np.array([u[pairs[:,0]], v[pairs[:,0]], w[pairs[:,0]]]).T
  u2 = np.array([u[pairs[:,1]], v[pairs[:,1]], w[pairs[:,1]]]).T
  u12 = u2 - u1

  # Sanity check using 3^2 + 4^2 + 12^2 = 13^2
  #u1[:,0] = 0; u1[:,1] = 0; u1[:,2] = 0;
  #u2[:,0] = 3; u2[:,1] = 4; u2[:,2] = 12;
  #u12 = u2 - u1

  # Find velocity difference perpendicular and parallel to:
  # -- vertical, ez = [0,0,1]
  # -- horizontal, exy = [1,1,0]/sqrt(2)
  du_verti_temp = u12[:,2]
  du_horiz_temp = linalg.norm(u12[:,(0,1)], axis=1)

  # -- separation vector, er = r12/r (parallel)
  du_parallel_temp = np.array([np.dot(uu,rvec/rr) for (uu, rvec, rr) in zip(u12, r12, r_mag)])

  # -- separation vector, u = u12 - dot(u_parallel, er)
  # -- separation vector, u = sqrt(u12**2 - u_parallel**2)
  du_perp_temp = np.sqrt(np.linalg.norm(u12, axis=1)**2 - du_parallel_temp**2)

  # Square everything for second order structure function
  du_verti_temp *= du_verti_temp
  du_horiz_temp *= du_horiz_temp
  du_parallel_temp *= du_parallel_temp
  du_perp_temp *= du_perp_temp

  # Sort the velocities by separation r
  index = np.argsort(r_mag)
  r_mag = r_mag[index]
  du_verti_temp = du_verti_temp[index]
  du_horiz_temp = du_horiz_temp[index]
  du_parallel_temp = du_parallel_temp[index]
  du_perp_temp = du_perp_temp[index]

  # Loop through bins and increment sum and counter
  for rr,center in enumerate(bin_center):
    # Find all indices that fall inside the bin 
    greater_ind = np.argwhere(r_mag >= (center - 0.5*bin_width))
    less_ind = np.argwhere(r_mag <= (center + 0.5*bin_width))
    ind = np.intersect1d(greater_ind, less_ind)

    # Increment counter
    du_bin_counts[rr] += len(ind)

    # Add sum
    du_verti[rr] += np.sum(du_verti_temp[ind])
    du_horiz[rr] += np.sum(du_horiz_temp[ind])
    du_parallel[rr] += np.sum(du_parallel_temp[ind])
    du_perpendicular[rr] += np.sum(du_perp_temp[ind])

  x0 = x
  y0 = y
  z0 = z

  bbparts.close()
print()

# Find mean
du_verti /= du_bin_counts
du_horiz /= du_bin_counts
du_parallel /= du_bin_counts
du_perpendicular /= du_bin_counts

## Temporary plotting ##
# Sort the velocities by r
#index = np.argsort(r_mag)
#r_mag = r_mag[index]
#u_z = u_z[index]
#u_xy = u_xy[index]
#u_sep_para = u_sep_para[index]
#u_sep_perp = u_sep_perp[index]

# Bin the separation distrance
#dr = 0.1*a           # Distance between bin centers
#bin_width = 0.2*a    # Width of bin. Should be >= dr, 
#bin_center = np.arange(2.*a + 0.5*dr, r, dr)
#nbins = len(bin_center)

#u_z_bin = np.zeros(nbins)
#u_xy_bin = np.zeros(nbins)
#u_para_bin = np.zeros(nbins)
#u_perp_bin = np.zeros(nbins)
#for rr,center in enumerate(bin_center):
#  # Find all indices that fall inside the bin 
#  greater_ind = np.argwhere(r_mag >= (center - 0.5*bin_width))
#  less_ind = np.argwhere(r_mag <= (center + 0.5*bin_width))
#  ind = np.intersect1d(greater_ind, less_ind)
#
#  u_z_bin[rr] = np.mean(u_z[ind])
#  u_xy_bin[rr] = np.mean(u_xy[ind])
#  u_para_bin[rr] = np.mean(u_sep_para[ind])
#  u_perp_bin[rr] = np.mean(u_sep_perp[ind])

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,6))

ax1 = fig.add_subplot(211)
plt.loglog(bin_center/a, du_verti, 'o', markersize=2, color='blue')
plt.loglog(bin_center/a, du_horiz, 'o', markersize=2, color='green')

plt.xlim([1.9, 10*r/a])
plt.ylim([1e-6, 5e-2])
plt.ylabel(r"Sg")
plt.legend([r'z', r'xy'])

ax1 = fig.add_subplot(212)
plt.loglog(bin_center/a, du_parallel, 'o', markersize=2, color='blue')
plt.loglog(bin_center/a, du_perpendicular, 'o', markersize=2, color='green')

plt.xlim([1.9, 10*r/a])
plt.ylim([5e-7, 1e-2])
plt.xlabel(r"r/a")
plt.ylabel(r"Ssep")
plt.legend([r'para', r'perp'])

plt.tight_layout()
plt.savefig("test-structure.png", bbox_inches='tight', format='png')
