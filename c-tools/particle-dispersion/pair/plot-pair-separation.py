#!/usr/bin/env python

# PURPOSE
#   Calculate the autocorrelation of particle velocity from the supplied time
#   series.
# USAGE
#   ./plot-pair-separation </path/to/sim/output>
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
  print("   ./plot-pair-separation.py <./path/to/sim/output>")
  sys.exit()

# Get data
savefile = data_dir + "/../analysis/pair-dispersion/data/pair-separation.csv"

data = np.genfromtxt(savefile, skip_header=1)

time = data[:,0]
r2_total     = data[:,1]
r2_verti     = data[:,2]
r2_horiz     = data[:,3]
cos_theta_z  = data[:,4]
u2_total     = data[:,5]
u2_verti     = data[:,6]
u2_horiz     = data[:,7]
u2_paral     = data[:,8]
u2_perpe     = data[:,9]

# Predefine simulation parameters -- XXX should not hardcode!
a = 2.1
Lx = 42.
Ly = 42.
Lz = 126.
nu = 0.01715  # [mm^2/ms]
rho = 3.3     # TODO pull from sim directory
tau_p = (2.*a)**2 * rho / (18. * nu)  # [ms]

# Set up imgidr
imgdir = data_dir + "/../analysis/pair-dispersion/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Separation #
plt.figure()
plt.loglog(time/tau_p, np.sqrt(r2_total)/a, '.', markersize=1)
plt.loglog(time/tau_p, np.sqrt(r2_verti)/a, '.', markersize=1)
plt.loglog(time/tau_p, np.sqrt(r2_horiz)/a, '.', markersize=1)

plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$\sqrt{\langle r_{\alpha \beta}^2\rangle}/a$")
plt.legend(["total", "verti", "horiz"])

plt.ylim(ymin=1)

# r ~ t^1/2
xpts = [1e0, 1e2]
ypts = np.power(xpts, 0.5)*10.
plt.loglog(xpts, ypts, 'k--')
plt.text(1e0, 20, r"$t^{0.5}$")

# r ~ t
#xpts = [2e-1, 1e0]
#ypts = np.power(xpts, 1)*10
#plt.loglog(xpts, ypts, 'k--')
#plt.text(1e0, 2, r"$t^1$")

plt.tight_layout()
plt.savefig(imgdir + "pair-separation.png", bbox_inches='tight', format='png')

# Alignment with gravity #
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.semilogx(time/tau_p, cos_theta_z, '.', markersize=1)
ax1.set_xlabel(r"$t/\tau_p$")
ax1.set_ylabel(r"$\langle \cos(\theta_z) \rangle_{n_p}$")
ax1.set_ylim([0, 1])

plt.tight_layout()
plt.savefig(imgdir + "pair-alignment.png", bbox_inches='tight', format='png')

# Velocity #
plt.figure(figsize=(5,2.5))
plt.plot(time/tau_p, np.sqrt(u2_verti), '-')#, markersize=1)
plt.plot(time/tau_p, np.sqrt(u2_horiz), '-')#, markersize=1)
plt.plot(time/tau_p, np.sqrt(u2_paral), '-')#, markersize=1)
plt.plot(time/tau_p, np.sqrt(u2_perpe), '-')#, markersize=1)
plt.plot(time/tau_p, np.sqrt(u2_total), '-')#, markersize=1)

plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$\sqrt{\langle u_{\alpha \beta}^2\rangle}$")
plt.legend(["verti", "horiz", 'parallel', 'perpendicular', "total"])

plt.xlim([0, 10])

plt.gca().set_prop_cycle(None)
plt.axhline(np.mean(np.sqrt(u2_total)), linestyle='--')
plt.axhline(np.mean(np.sqrt(u2_verti)), linestyle='--')
plt.axhline(np.mean(np.sqrt(u2_horiz)), linestyle='--')
plt.axhline(np.mean(np.sqrt(u2_paral)), linestyle='--')
plt.axhline(np.mean(np.sqrt(u2_perpe)), linestyle='--')

plt.tight_layout()
plt.savefig(imgdir + "pair-velocity.png", bbox_inches='tight', format='png')


