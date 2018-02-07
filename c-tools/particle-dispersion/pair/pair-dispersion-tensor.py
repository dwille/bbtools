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

r = 3.*a    # kd-tree search radius


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
npairs = len(pairs)
print("  Found %d pairs within a distance of %.2fa" % (npairs, r/a))

### Init memory ###
flip_i = np.zeros(nparts)
flip_j = np.zeros(nparts)
flip_k = np.zeros(nparts)

D = np.zeros((3,3,len(times)))

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

  ### Dispersion Tensor ###
  # From Shen and Yeung
  # D_ij(t) = < (r_i^(1) - r_i^(2))(r_j^(1) - r_j^(2)) >

  # Find pair separtion. r1 is position of 1st particle, r2 is position of 2nd
  #r1 = np.array([x[pairs[:,0]], y[pairs[:,0]], z[pairs[:,0]]]).T
  #r2 = np.array([x[pairs[:,1]], y[pairs[:,1]], z[pairs[:,1]]]).T
  #r12 = r1 - r2
  # Vectorize this?

  # Find dispersion tensor for each particle
  for pp in np.arange(npairs):
    p1 = pairs[pp,0]
    p2 = pairs[pp,1]

    r1 = np.matrix([x[p1], y[p1], z[p1]])
    r2 = np.matrix([x[p2], y[p2], z[p2]])
    r12 = r1 - r2

    # r12.T * r12 multiples a (3x1)*(1x3) to (3x3)
    # D[:,:,tt] = r12.T * r12
    D[:,:,tt] += r12.T * r12
  
  D[:,:,tt] /= npairs

  # Save data as old timestep
  x0 = x
  y0 = y
  z0 = z

  bbparts.close()
print()

# Save data to file
import csv,os
savefile = data_dir + "/../analysis/pair-dispersion/data/"
if not os.path.exists(savefile):
  os.makedirs(savefile)
savefile += "pair-dispersion-tensor.csv"

with open(savefile, 'w') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerow(["time","D11","D12","D13","D21","D22","D23","D31","D32","D33"])
  data_out = np.zeros((len(time_flt),10))
  data_out[:,0] = time_flt
  data_out[:,1] = D[0,0,:]
  data_out[:,2] = D[0,1,:]
  data_out[:,3] = D[0,2,:]
  data_out[:,4] = D[1,0,:]
  data_out[:,5] = D[1,1,:]
  data_out[:,6] = D[1,2,:]
  data_out[:,7] = D[2,0,:]
  data_out[:,8] = D[2,1,:]
  data_out[:,9] = D[2,2,:]
  out.writerows(data_out)
