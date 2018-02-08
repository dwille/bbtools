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
import bluebottle_particle_reader as bbparts
from matplotlib.ticker import MultipleLocator

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

# Set up imgidr
imgdir = data_dir + "/../analysis/pair-dispersion/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

##################################
########### Separation ###########
##################################

plt.figure()
plt.loglog(time/tau_p, np.sqrt(r2_total)/a, '.', markersize=1)
plt.loglog(time/tau_p, np.sqrt(r2_verti)/a, '.', markersize=1)
plt.loglog(time/tau_p, np.sqrt(r2_horiz)/a, '.', markersize=1)

plt.xlabel(r"$t/\tau_p$")
plt.ylabel(r"$\sqrt{\langle r_{\alpha \beta}^2\rangle}/a$")
plt.legend(["total", "verti", "horiz"])
plt.title(r"$\rho^* = %.1f,\ N_p = %d$" % (rho, nparts))

plt.xlim([1e-2,2e2])
plt.ylim([1,100])

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

# mean free path
lambda_1 = 1. / (3. * phi)
lambda_2 = 4. / (27. * phi)

plt.axhline(np.sqrt(r2_total[0])/a + lambda_1, linestyle='--', color='k')
plt.axhline(np.sqrt(r2_total[0])/a + lambda_2, linestyle='--', color='k')

plt.tight_layout()
plt.savefig(imgdir + "pair-separation.png", bbox_inches='tight', format='png')


##############################################
########### Alignment with gravity ###########
##############################################

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.semilogx(time/tau_p, cos_theta_z, '.', markersize=1)
ax1.set_xlabel(r"$t/\tau_p$")
ax1.set_ylabel(r"$\langle \cos(\theta_z) \rangle_{n_p}$")
ax1.set_ylim([0, 1])

plt.title(r"$\rho^* = %.1f,\ N_p = %d$" % (rho, nparts))

plt.tight_layout()
plt.savefig(imgdir + "pair-alignment.png", bbox_inches='tight', format='png')

################################
########### Velocity ###########
################################

fig2 = plt.figure(figsize=(4,4))

# Aligned with csys
ax2a = fig2.add_subplot(211)
ax2a.plot(time/tau_p, np.sqrt(u2_verti), '-')
ax2a.plot(time/tau_p, np.sqrt(u2_horiz), '-')
ax2a.plot(time/tau_p, np.sqrt(u2_total), '-')

ax2a.set_xlabel(r"$t/\tau_p$")
ax2a.set_ylabel(r"$\sqrt{\langle u_{\alpha \beta}^2\rangle}$")
plt.title(r"$\rho^* = %.1f,\ N_p = %d$" % (rho, nparts))

ax2a.legend([r"$u_z$", r"$u_{xy}$", r"$\left\Vert \mathbf{u} \right\Vert$"],
            loc='right')

ax2a.set_xlim([0, 10])
ax2a.set_ylim(ymin=0)
ax2a.xaxis.set_major_locator(MultipleLocator(2))
ax2a.xaxis.set_minor_locator(MultipleLocator(0.5))

# Aligned with spearation
ax2b = fig2.add_subplot(212)
ax2b.plot(time/tau_p, np.sqrt(u2_paral), '-')
ax2b.plot(time/tau_p, np.sqrt(u2_perpe), '-')
ax2b.plot(time/tau_p, np.sqrt(u2_total), '-')

ax2b.set_xlabel(r"$t/\tau_p$")
ax2b.set_ylabel(r"$\sqrt{\langle u_{\alpha \beta}^2\rangle}$")

ax2b.legend([r"$u_\parallel$", r"$u_\perp$", r"$\left\Vert \mathbf{u} \right\Vert$"],
             loc='right')

ax2b.set_xlim([0, 10])
ax2b.set_ylim(ymin=0)
ax2b.xaxis.set_major_locator(MultipleLocator(2))
ax2b.xaxis.set_minor_locator(MultipleLocator(0.5))

# need to import some other modules to cycle color correctly here, since 
# axhline doesn't seem to by default
# https://stackoverflow.com/questions/30535442/matplotlib-fill-between-does-not-cycle-through-colours
# https://stackoverflow.com/questions/34247297/matplotlib-1-5-usage-of-axes-prop-cycle#34248178
from matplotlib import rcParams
import itertools

ax2a.set_prop_cycle(None)
clist = rcParams['axes.color_cycle']
cgen = itertools.cycle(clist)

ax2a.axhline(np.mean(np.sqrt(u2_verti)), linestyle='--', color=next(cgen))
ax2a.axhline(np.mean(np.sqrt(u2_horiz)), linestyle='--', color=next(cgen))
ax2a.axhline(np.mean(np.sqrt(u2_total)), linestyle='--', color=next(cgen))

ax2b.set_prop_cycle(None)
clist = rcParams['axes.color_cycle']
cgen = itertools.cycle(clist)

ax2b.axhline(np.mean(np.sqrt(u2_paral)), linestyle='--', color=next(cgen))
ax2b.axhline(np.mean(np.sqrt(u2_perpe)), linestyle='--', color=next(cgen))
ax2b.axhline(np.mean(np.sqrt(u2_total)), linestyle='--', color=next(cgen))

plt.tight_layout()
plt.savefig(imgdir + "pair-velocity.png", bbox_inches='tight', format='png')


