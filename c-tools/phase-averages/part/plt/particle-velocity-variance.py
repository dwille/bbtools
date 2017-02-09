#!/usr/bin/env python2
from setup import *

os.system('clear')
print ""
print " ---- Phase-Averaged Particle Velocity Plotting Utility ---- "
print ""

# Parse command line args and set up directory structure
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/"
imgdir = root + "simdata/img/part/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)
  print "Made directory %s" % imgdir

# Sim params XXX
partR = 2.1
domX = 42
domY = 42
domZ = 126
nparts = np.array([500, 1000, 1500, 2000])
rho = np.array([2.0, 3.3, 4.0, 5.0])
phi = 4./3.*np.pi*partR*partR*partR*nparts / (domX*domY*domZ)

# Files
varFile = root + "simdata/part_vel_variance"
wfFile = root + "simdata/phaseAveragedFluidVel"
# Check if terminal vel file exists
if not os.path.exists(varFile):
  print "      " + varFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(wfFile):
  print "      " + wfFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Pull data
vert_var = np.genfromtxt(varFile, skip_header=1, usecols=2)
hori_var = np.genfromtxt(varFile, skip_header=1, usecols=3)
var_ratio = vert_var / hori_var

phase_avg_wf = np.genfromtxt(wfFile, skip_header=1, usecols=2)
vert_var /= phase_avg_wf
hori_var /= phase_avg_wf

# Plot
fig1 = plt.figure()
ax1 = fig1.add_subplot(311)

ax1.plot(phi, var_ratio[0:16:4], 'ko--')
ax1.plot(phi, var_ratio[1:16:4], 'ro--')
ax1.plot(phi[0:3], var_ratio[2:12:4], 'go--')
ax1.plot(phi[0:2], var_ratio[3:8:4], 'bo--')

ax1.set_xlim([0, 0.4])
ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.set_xticklabels([])

ax1.set_ylim([0, 2.5])
ax1.set_ylabel(r"$\frac{\mathrm{Var}[u_z]}{\mathrm{Var}[u_\perp]}$",
  rotation=0)
ax1.yaxis.set_label_coords(-0.2, 0.5)
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(MultipleLocator(0.25))

ax1.legend([r"$\rho^* = 2$",r"$\rho^* = 3.3$",r"$\rho^* = 4$",r"$\rho^* = 5$"],
  loc="lower left", fontsize=6, framealpha=0.7)

#
ax2 = fig1.add_subplot(312)

ax2.plot(phi, vert_var[0:16:4], 'ko--')
ax2.plot(phi, vert_var[1:16:4], 'ro--')
ax2.plot(phi[0:3], vert_var[2:12:4], 'go--')
ax2.plot(phi[0:2], vert_var[3:8:4], 'bo--')

ax2.plot(phi, hori_var[0:16:4], 'ks-.')
ax2.plot(phi, hori_var[1:16:4], 'rs-.')
ax2.plot(phi[0:3], hori_var[2:12:4], 'gs-.')
ax2.plot(phi[0:2], hori_var[3:8:4], 'bs-.')

ax2.set_xlim([0, 0.4])
ax2.xaxis.set_major_locator(MultipleLocator(0.1))
ax2.xaxis.set_minor_locator(MultipleLocator(0.05))
ax2.set_xticklabels([])

ax2.set_ylim([0, .5])
ax2.set_ylabel(r"$\frac{\mathrm{Var}[u_i]}{\langle w_f \rangle}$",
  rotation=0)
ax2.yaxis.set_label_coords(-0.2, 0.5)
ax2.yaxis.set_major_locator(MultipleLocator(0.25))
ax2.yaxis.set_minor_locator(MultipleLocator(0.125))

#
ax3 = fig1.add_subplot(313)

ax3.plot(phi, phase_avg_wf[0:16:4], 'ko--')
ax3.plot(phi, phase_avg_wf[1:16:4], 'ro--')
ax3.plot(phi, phase_avg_wf[2:16:4], 'go--')
ax3.plot(phi, phase_avg_wf[3:16:4], 'bo--')

ax3.set_xlim([0, 0.4])
ax3.set_xlabel(r"$\phi$")
ax3.xaxis.set_major_locator(MultipleLocator(0.1))
ax3.xaxis.set_minor_locator(MultipleLocator(0.05))

ax3.set_ylim([0, .5])
ax3.set_ylabel(r"$\langle w_f \rangle$", rotation=0)
ax3.yaxis.set_label_coords(-0.2, 0.5)
ax3.yaxis.set_major_locator(MultipleLocator(0.25))
ax3.yaxis.set_minor_locator(MultipleLocator(0.125))



imgname = imgdir + "vel_variance"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
