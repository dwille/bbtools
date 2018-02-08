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
import matplotlib.pyplot as plt
import bluebottle_particle_reader as bbparts
from matplotlib.ticker import MultipleLocator

##########

# Parse output directory from command line
if len(sys.argv) >= 1:    # output directory given
  data_dir = sys.argv[1]

else:                     # nothing given
  print("Error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./plot-part-pair-separation.py <./path/to/sim/output>")
  sys.exit()

# Get data
savefile = data_dir + "/../analysis/pair-dispersion/data/pair-vel-autocorrelation.csv"

data = np.genfromtxt(savefile, skip_header=1)

time = data[:,0]
u_verti        = data[:,1]
u_verti_horiz  = data[:,2]
u_horiz        = data[:,3]
u_paral        = data[:,4]
u_paral_perpe  = data[:,5]
u_perpe        = data[:,6]

# Normalize to one
#u_verti /= u_verti[0]
#u_horiz /= u_horiz[0]
#u_paral /= u_paral[0]
#u_perpe /= u_perpe[0]

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

# Find response time based on d^2R/dt^2 at origin
# d2Rdt2 = R_{i+1} + R_{i-1} - 2R_i / dt^2
# but R is even over t = 0, so R_{i+-1} are equal
# d2Rdt2 = 2R_{i+1} - 2R_i / dt^2 -> b < 0
# dRdt = bt + c1
# dRdt(t=0) = 0 -> c1 =0
# R = 0.5*bt^2 + c2
# R(t=0) = c2 = R0
# -> R = 0.5*bt^2 + R0
# Want t where R = 0
#   t = sqrt( -2*R0 / b )

dt = (time[1] - time[0])/tau_p

diff2 = 2.*(u_verti[1] - u_verti[0]) / dt**2
t_verti = np.sqrt(-2. * u_verti[0] / diff2)

diff2 = 2.*(u_horiz[1] - u_horiz[0]) / dt**2
t_horiz = np.sqrt(-2. * u_horiz[0] / diff2)

diff2 = 2.*(u_paral[1] - u_paral[0]) / dt**2
t_paral = np.sqrt(-2. * u_paral[0] / diff2)

diff2 = 2.*(u_perpe[1] - u_perpe[0]) / dt**2
t_perpe = np.sqrt(-2. * u_perpe[0] / diff2)

# Set up imgidr
imgdir = data_dir + "/../analysis/pair-dispersion/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Plot
fig = plt.figure(figsize=(4,4))

# Aligned with coord sys
ax1 = fig.add_subplot(311)
plt.plot(time/tau_p, u_verti)
plt.plot(time/tau_p, u_horiz)

# Use parabola to find estimate of correlation time
diff2 = 2.*(u_verti[1] - u_verti[0]) / dt**2
tplot = np.arange(0, t_verti, dt/5.)
plt.plot(tplot , 0.5*diff2*tplot**2 + u_verti[0], 'k--', linewidth=1)

diff2 = 2.*(u_horiz[1] - u_horiz[0]) / dt**2
tplot = np.arange(0, t_horiz, dt/5.)
plt.plot(tplot , 0.5*diff2*tplot**2 + u_horiz[0], 'k--', linewidth=1)

plt.legend([r"$u_z$", r"$u_{xy}$"])
plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$R_{ff}(\tau)$")
plt.title(r"$\rho^* = %.1f,\ N_p = %d$" % (rho, nparts))
plt.xlim([0, 2])
plt.ylim(ymin=0)

ax1.xaxis.set_major_locator(MultipleLocator(0.5))
ax1.xaxis.set_minor_locator(MultipleLocator(0.1))

ax2 = fig.add_subplot(312)
# Aligned with separation
plt.plot(time/tau_p, u_paral)
plt.plot(time/tau_p, u_perpe)

# Use parabola to find estimate of correlation time
diff2 = 2.*(u_paral[1] - u_paral[0]) / dt**2
tplot = np.arange(0, t_paral, dt/5.)
plt.plot(tplot , 0.5*diff2*tplot**2 + u_paral[0], 'k--', linewidth=1)

diff2 = 2.*(u_perpe[1] - u_perpe[0]) / dt**2
tplot = np.arange(0, t_perpe, dt/5.)
plt.plot(tplot , 0.5*diff2*tplot**2 + u_perpe[0], 'k--', linewidth=1)

plt.legend([r"$u_\parallel$", r"$u_\perp$"])
plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$R_{ff}(\tau)$")
plt.xlim([0, 2])
plt.ylim(ymin=0)

ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.xaxis.set_minor_locator(MultipleLocator(0.1))

ax3 = fig.add_subplot(313)
# Cross-correlations
plt.plot(time/tau_p, u_verti_horiz)
plt.plot(time/tau_p, u_paral_perpe)

plt.legend([r"$u_{z,xy}$", r"$u_{\parallel,\perp}$"],
            loc='upper right')
plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$R_{fg}(\tau)$")
plt.xlim([0, 2])
plt.ylim(ymin=0)

ax3.xaxis.set_major_locator(MultipleLocator(0.5))
ax3.xaxis.set_minor_locator(MultipleLocator(0.1))

plt.tight_layout()
plt.savefig(imgdir + "pair-vel-autocorrelation.png", bbox_inches='tight', format='png')

