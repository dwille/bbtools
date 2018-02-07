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


# Predefine simulation parameters -- XXX should not hardcode!
a = 2.1
Lx = 42.
Ly = 42.
Lz = 126.
nu = 0.01715  # [mm^2/ms]
rho = 3.3     # TODO pull from sim directory
tau_p = (2.*a)**2 * rho / (18. * nu)  # [ms]

# Plot

# Set up imgidr
imgdir = data_dir + "/../analysis/pair-dispersion/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

fig = plt.figure()

ax1 = fig.add_subplot(211)
plt.plot(time/tau_p, u_verti) #/u_verti[0])
plt.plot(time/tau_p, u_horiz) #/u_horiz[0])
plt.plot(time/tau_p, u_paral) #/u_paral[0])
plt.plot(time/tau_p, u_perpe) #/u_perpe[0])

plt.legend([r"$u_z$", r"$u_{xy}$", r"$u_\parallel$", r"$u_\perp$"])
plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$R_{ff}(\tau)$")
plt.xlim([0, 2])

ax2 = fig.add_subplot(212)
plt.plot(time/tau_p, u_verti_horiz)
plt.plot(time/tau_p, u_paral_perpe)

plt.legend([r"$u_{z,xy}$", r"$u_{\parallel,\perp}$"])
plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$R_{fg}(\tau)$")
plt.xlim([0, 2])

plt.tight_layout()
plt.savefig(imgdir + "pair-vel-autocorrelation.png", bbox_inches='tight', format='png')

