#!/usr/bin/env python

# PURPOSE
#   Calculate the autocorrelation of particle velocity from the supplied time
#   series.
# USAGE
#   ./part-pair-separation.py </path/to/sim/output> <start_time>
#
# OUTPUT
#   </path/to/sim/analysis/pair-dispersion/img>

# Imports:
import sys
import numpy as np
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
nu = 0.01715  # [mm^2/ms]
rho = 3.3     # TODO pull from sim directory
tau_p = (2.*a)**2 * rho / (18. * nu)  # [ms]

r = 3.*a      # kd-tree search radius

####################################################
#### Initialize particle pairs using a k-d-tree ####
####################################################

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

# Set up tree with particle positions
data = np.c_[x0.ravel(), y0.ravel(), z0.ravel()]
tree = spatial.cKDTree(data)#, boxsize=boxsize)

# Query particle pairs less than desired radius
pairs = tree.query_pairs(r)
pairs = np.array(list(pairs))
p1 = pairs[:,0]
p2 = pairs[:,1]
npairs = len(pairs)
print("  Found %d pairs within a distance of %.2f" % (npairs, r))


# Init memory 
r2_total = np.zeros(len(times))
r2_horiz = np.zeros(len(times))
r2_verti = np.zeros(len(times))
cos_theta_z = np.zeros(len(times))

u2_total = np.zeros(len(times))
u2_horiz = np.zeros(len(times))
u2_verti = np.zeros(len(times))
u2_paral = np.zeros(len(times))
u2_perpe = np.zeros(len(times))

flip_i = np.zeros(nparts)
flip_j = np.zeros(nparts)
flip_k = np.zeros(nparts)

### Loop over time and find mean squared separation ###
for tt,time in enumerate(times):
  print("  Time %d of %d\r" % (tt+1, len(times)), end="")
  bbparts.open(time)
  (x, y, z) = bbparts.read_part_position();
  (u, v, w) = bbparts.read_part_velocity();

  ### Deal with particles crossing periodic boundaries ###
  if (tt > 0):
   (x, y, z, flip_i, flip_j, flip_k) =                      \
                bbparts.periodic_flip(x, y, z, x0, y0, z0,  \
                             flip_i, flip_j, flip_k,        \
                             Lx, Ly, Lz)

  ## Particle psotions ###
  r1 = np.array([x[p1], y[p1], z[p1]]).T
  r2 = np.array([x[p2], y[p2], z[p2]]).T
  r12 = r1 - r2

  ## Find mean-squared separation
  squared = np.power(np.linalg.norm(r12, axis=1), 2)
  r2_total[tt] = np.mean(squared)

  squared = np.power(r12[:,2], 2)
  r2_verti[tt] = np.mean(squared)

  squared = np.power(np.linalg.norm(r12[:,(0,1)], axis=1), 2)
  r2_horiz[tt] = np.mean(squared)

  # Alignment with vertical
  # cos(theta_z) = |r_z| / ||r12||
  # Not concerned with up/down, so do |dz|
  arg = np.abs(r12[:,2]) / np.linalg.norm(r12, axis=1)
  cos_theta_z[tt] = np.mean(arg)

  ## Particle Velocities ##
  u1 = np.array([u[p1], v[p1], w[p1]]).T
  u2 = np.array([u[p2], v[p2], w[p2]]).T
  u12 = u1 - u2

  #...xyz coordinates
  squared = np.power(np.linalg.norm(u12, axis=1), 2)
  u2_total[tt] = np.mean(squared)

  squared = np.power(u12[:,2], 2)
  u2_verti[tt] = np.mean(squared)

  squared = np.power(np.linalg.norm(u12[:,(0,1)], axis=1), 2)
  u2_horiz[tt] = np.mean(squared)

  # ...relative to separation
  e_sep = r12 / np.linalg.norm(r12, axis=1)[:,np.newaxis]

  u_parallel = np.sum(u12 * e_sep, axis=1)
  u2_paral[tt] = np.mean(u_parallel**2)

  u_squared = np.power(np.linalg.norm(u12, axis=1), 2)
  u2_perpe[tt] = np.mean(u_squared - u_parallel**2)

  # Save data as old timestep
  x0 = x
  y0 = y
  z0 = z

  bbparts.close()
print()

# Save data
import csv, os
savefile = data_dir + "/../analysis/pair-dispersion/data/"
if not os.path.exists(savefile):
  os.makedirs(savefile)
savefile += "pair-separation.csv"

with open(savefile, 'w') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerow(["time","r2_total", "r2_verti", "r2_horiz", "cos_theta", "u2_total", \
                 "u2_verti", "u2_horiz", "u2_paral", "u2_perpe"])
  data_out = np.zeros((len(times), 10))
  data_out[:,0] = time_flt
  data_out[:,1] = r2_total
  data_out[:,2] = r2_verti
  data_out[:,3] = r2_horiz
  data_out[:,4] = cos_theta_z
  data_out[:,5] = u2_total
  data_out[:,6] = u2_verti
  data_out[:,7] = u2_horiz
  data_out[:,8] = u2_paral
  data_out[:,9] = u2_perpe
  out.writerows(data_out)

