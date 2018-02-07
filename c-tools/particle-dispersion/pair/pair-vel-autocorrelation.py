#!/usr/bin/env python

# PURPOSE
#   Calculate the autocorrelation of particle velocity from the supplied time
#   series.
# Two-point correlation function: R_ij = < u_i(x + r, t) u_j(x, t)>
# USAGE
#   ./part-pair-separation.py </path/to/sim/output> <start_time>
#
# OUTPUT
#   </path/to/sim/analysis/pair-dispersion/img>

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
u_verti = np.zeros((len(times), npairs))
u_horiz = np.zeros((len(times), npairs))
u_paral = np.zeros((len(times), npairs))
u_perpe = np.zeros((len(times), npairs))

flip_i = np.zeros(nparts)
flip_j = np.zeros(nparts)
flip_k = np.zeros(nparts)

### Loop over time and pull data ###
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

  ## Particle separation ##
  r1 = np.array([x[p1], y[p1], z[p1]]).T
  r2 = np.array([x[p2], y[p2], z[p2]]).T
  r12 = r1 - r2

  e_sep = r12 / np.linalg.norm(r12, axis=1)[:,np.newaxis]

  ## Relative velocity
  u1 = np.array([u[p1], v[p1], w[p1]]).T
  u2 = np.array([u[p2], v[p2], w[p2]]).T
  u12 = u1 - u2

  ## Save
  u_verti[tt,:] = u12[:,2]
  u_horiz[tt,:] = np.sqrt(u12[:,0]**2 + u12[:,1]**2)

  u_paral[tt,:] = np.sum(u12 * e_sep, axis=1)

  u_squared = np.power(np.linalg.norm(u12, axis=1), 2)
  u_perpe[tt,:] = np.sqrt(u_squared - u_paral[tt,:]**2)

  # Save data as old timestep
  x0 = x
  y0 = y
  z0 = z

  bbparts.close()
print()

#### Perform two-point correlation function ####
# rho_ij = < u_i(t) u_j(t + tau)>
# have:
# 1) verti-verti 2) verti-horiz
# 3) horiz-horiz 4) paral-paral
# 5) paral-perpe 5) perpe-perpe
print("Performing autocorrelation...")

u_verti_horiz = np.zeros((len(times), npairs))
u_paral_perpe = np.zeros((len(times), npairs))
for nn in np.arange(npairs):
  print(" pair %d of %d\r" % (nn+1, npairs), end="")

  # Subtract mean from velocities
  t_u_verti = u_verti[:,nn] - np.mean(u_verti[:,nn])
  t_u_horiz = u_horiz[:,nn] - np.mean(u_horiz[:,nn])
  t_u_paral = u_paral[:,nn] - np.mean(u_paral[:,nn])
  t_u_perpe = u_perpe[:,nn] - np.mean(u_perpe[:,nn])

  # Convolutions
  # place autocorrelations back into the array
  # place crosscorrelations into another array

  # XXX why do we do [::-1]?

  # 2) verti-horiz
  result = signal.fftconvolve(t_u_verti[::-1], t_u_horiz, mode="full")
  result = result[::-1]
  pstart = int(np.floor(len(result)/2.))
  u_verti_horiz[:,nn] = result[pstart:]
  # normalize by setting R(tau=0)=1
  #u_verti_horiz[:,nn] /= u_verti_horiz[0,nn]
  # normalize by sdev(a)*sdev(b), sdev is over time
  #u_verti_horiz[:,nn] /= (np.std(t_u_verti)*np.std(t_u_horiz))

  # 1) verti-verti
  result = signal.fftconvolve(t_u_verti[::-1], t_u_verti, mode="full")
  result = result[::-1]
  pstart = int(np.floor(len(result)/2.))
  u_verti[:,nn] = result[pstart:]
  #u_verti[:,nn] /= u_verti[0,nn]
  # normalize by sdev(a)*sdev(b), sdev is over time
  #u_verti[:,nn] /= np.std(t_u_verti)**2

  # 3) horiz-horiz
  result = signal.fftconvolve(t_u_horiz[::-1], t_u_horiz, mode="full")
  result = result[::-1]
  pstart = int(np.floor(len(result)/2.))
  u_horiz[:,nn] = result[pstart:]
  #u_horiz[:,nn] /= u_horiz[0,nn]
  # normalize by sdev(a)*sdev(b), sdev is over time
  #u_horiz[:,nn] /= np.std(t_u_horiz)**2

  # 5) paral-perpe
  result = signal.fftconvolve(t_u_paral[::-1], t_u_perpe, mode="full")
  result = result[::-1]
  pstart = int(np.floor(len(result)/2.))
  u_paral_perpe[:,nn] = result[pstart:]
  #u_paral_perpe[:,nn] /= u_paral_perpe[0,nn]
  #u_paral_perpe[:,nn] /= (np.std(t_u_paral)*np.std(t_u_perpe))

  # 4) paral-paral
  result = signal.fftconvolve(t_u_paral[::-1], t_u_paral, mode="full")
  result = result[::-1]
  pstart = int(np.floor(len(result)/2.))
  u_paral[:,nn] = result[pstart:]
  #u_paral[:,nn] /= u_paral[0,nn]
  #u_paral[:,nn] /= np.std(t_u_paral)**2

  # 3) perpe-perpe
  result = signal.fftconvolve(t_u_perpe[::-1], t_u_perpe, mode="full")
  result = result[::-1]
  pstart = int(np.floor(len(result)/2.))
  u_perpe[:,nn] = result[pstart:]
  #u_perpe[:,nn] /= u_perpe[0,nn]
  #u_perpe[:,nn] /= np.std(t_u_perpe)**2

print() # for pretty output

# Average over each pair
u_verti = np.mean(u_verti, axis=1)
u_horiz = np.mean(u_horiz, axis=1)
u_verti_horiz = np.mean(u_verti_horiz, axis=1)

u_paral = np.mean(u_paral, axis=1)
u_perpe = np.mean(u_perpe, axis=1)
u_paral_perpe = np.mean(u_paral_perpe, axis=1)

# Save data to file
import csv,os
savefile = data_dir + "/../analysis/pair-dispersion/data/"
if not os.path.exists(savefile):
  os.makedirs(savefile)
savefile += "pair-vel-autocorrelation.csv"

with open(savefile, 'w') as outfile:
  out = csv.writer(outfile, delimiter=' ')
  out.writerow(["time","u_verti","u_verti_horiz","u_horiz","u_paral","u_paral_perpe","u_paral"])
  data_out = np.zeros((len(time_flt),7))
  data_out[:,0] = time_flt
  data_out[:,1] = u_verti
  data_out[:,2] = u_verti_horiz
  data_out[:,3] = u_horiz
  data_out[:,4] = u_paral
  data_out[:,5] = u_paral_perpe
  data_out[:,6] = u_perpe
  out.writerows(data_out)


# temporary plot
#import matplotlib.pyplot as plt
