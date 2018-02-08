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
import sys,os
import numpy as np
import bluebottle_particle_reader as bbparts
import matplotlib.pyplot as plt

##########

# Parse output directory from command line
if len(sys.argv) >= 1:    # output directory given
  data_dir = sys.argv[1]

else:                     # nothing given
  print("Error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./plot-part-pair-dispersion-tensor.py <./path/to/sim/output>")
  sys.exit()

# Get data
savefile = data_dir + "/../analysis/pair-dispersion/data/pair-dispersion-tensor.csv"

data = np.genfromtxt(savefile, skip_header=1)

time = data[:,0]
D = np.zeros((3,3,len(time)))
D[0,0,:] = data[:,1]
D[0,1,:] = data[:,2]
D[0,2,:] = data[:,3]
D[1,0,:] = data[:,4]
D[1,1,:] = data[:,5]
D[1,2,:] = data[:,6]
D[2,0,:] = data[:,7]
D[2,1,:] = data[:,8]
D[2,2,:] = data[:,9]

# Predefine simulation parameters
a = 2.1         # [mm]
Lx = 42.        # [mm]
Ly = 42.        # [mm]
Lz = 126.       # [mm]
nu = 0.01715    # [mm^2/ms]
rho_f = 8.75e-4 # [g/mm^3]

# Pull particle density from cgns file
cgns_times = bbparts.init(data_dir)
bbparts.open(cgns_times[0])
rho_p = bbparts.read_mean_part_density()
nparts = bbparts.read_nparts()

# Derived quantities
part_vol = 4./3. * np.pi * a**3
phi = nparts * part_vol / (Lx * Ly * Lz)
rho = rho_p / rho_f
tau_p = (2.*a)**2 * rho / (18. * nu)  # [ms]

# Plot

# Set up imgidr
imgdir = data_dir + "/../analysis/pair-dispersion/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

plt.figure()

plt.loglog(time / tau_p, np.sqrt(D[0,0,:])/a)
plt.loglog(time / tau_p, np.sqrt(D[1,1,:])/a)
plt.loglog(time / tau_p, np.sqrt(D[2,2,:])/a)

# ~ t^(1/2)
xpts = [1e1, 1e2]
ypts = np.power(xpts, 0.5)*2
plt.loglog(xpts, ypts, 'k--')
plt.text(xpts[0], 0.75*ypts[0], r"$t^{1/2}$")

# ~ t^1
xpts = [3e-1, 1e0]
ypts = np.power(xpts, 1.0)*10
plt.loglog(xpts, ypts, 'k--')
plt.text(xpts[0], 1.25*ypts[0], r"$t^1$")

# ~ t^(2/3)
xpts = [1.5e0, 2e1]
ypts = np.power(xpts, 2./3.)*10
plt.loglog(xpts, ypts, 'k--')
plt.text(xpts[0], 1.25*ypts[0], r"$t^{2/3}$")

plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$\sqrt{D_{ii}}/a$")
plt.ylim([1, 100])
plt.legend([r"$D_{11}$", r"$D_{22}$", r"$D_{33}$"])

plt.title("rho %.2f, nparts %d" % (rho, nparts))

plt.tight_layout()
plt.savefig(imgdir + "pair-dispersion-tensor.png", bbox_inches='tight', format='png')
