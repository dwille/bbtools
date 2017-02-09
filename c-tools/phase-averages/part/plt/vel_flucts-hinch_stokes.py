#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os,sys
os.system('clear')


print ""
print " ---- Phase-Averaged Particle Velocity Flucts Plotting Utility"
print ""

# Simulation Params
d = 4.2
nu = .01715
rhof = 8.75e-4
partR = 2.1
g = 0.00981
rhopf = np.array([2.0, 3.3, 4.0, 5.0])
rhop = np.array([0.001750, .00289, 0.0035, 0.004375])
phi = np.array([0.087, 0.175, 0.262, 0.349])

# MARCC
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
fluctFile = root + "velFlucts"
termFile = root + "singlePartSedi"

print "      Root directory set to: " + root
print "      Using location: " + fluctFile
print "      Using location: " + termFile
print ""
print "      Using d = %.2f" % d
print "      Using nu = %.5f" % nu

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(fluctFile):
  print "      " + fluctFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

imgdir = imgdir + "velocity_flucts/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Read data
uflucts_wrel = np.genfromtxt(fluctFile, skip_header=1, usecols=2)
vflucts_wrel = np.genfromtxt(fluctFile, skip_header=1, usecols=3)
wflucts_wrel = np.genfromtxt(fluctFile, skip_header=1, usecols=4)
wrel = np.genfromtxt(fluctFile)
dwdphi = np.genfromtxt(fluctFile, skip_header=1, usecols=8)

termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)
stokes_vel = (2.*rhopf*partR*partR*g)/(9.*nu)

normVel = stokes_vel

# Constant rho -- phi increasing
rho2_0u = uflucts_wrel[0:16:4]
rho2_0v = vflucts_wrel[0:16:4]
rho2_0w = wflucts_wrel[0:16:4]
rho3_3u = uflucts_wrel[1:16:4]
rho3_3v = vflucts_wrel[1:16:4]
rho3_3w = wflucts_wrel[1:16:4]
rho4_0u = uflucts_wrel[2:16:4]
rho4_0v = vflucts_wrel[2:16:4]
rho4_0w = wflucts_wrel[2:16:4]
rho5_0u = uflucts_wrel[3:16:4]
rho5_0v = vflucts_wrel[3:16:4]
rho5_0w = wflucts_wrel[3:16:4]

dwdphi2_0 = dwdphi[0:16:4]
dwdphi3_3 = dwdphi[1:16:4]
dwdphi4_0 = dwdphi[2:16:4]
dwdphi5_0 = dwdphi[3:16:4]

# Theoretical relations phi^1/3
# Hinch, 1988, Sedimentation of small particles (from disorder and mixing)
phiEval = np.linspace(0.01, 0.45, 81)
#highRe_20 = (normVel[0]*partR/nu)**(-1./3.)
#highRe_33 = (normVel[1]*partR/nu)**(-1./3.)
#highRe_40 = (normVel[2]*partR/nu)**(-1./3.)
#highRe_50 = (normVel[3]*partR/nu)**(-1./3.)
const = 0.5

## Plot Constant rho -- phi increasing ##
fig1 = plt.figure(figsize=(4,3))

ax1 = fig1.add_subplot(211)
ax1.plot(phi, rho2_0w, 'b*', markersize=7, alpha=0.7)
ax1.plot(phi, rho3_3w, 'gs', markersize=7, alpha=0.7)
ax1.plot(phi, rho4_0w, 'ro', markersize=7, alpha=0.7)
ax1.plot(phi, rho5_0w, 'c^', markersize=7, alpha=0.7)

#ax1.set_ylabel(r'$w^* = \sigma_{w_p}/\bar{w}_{rel}$', fontsize=14)
ax1.set_ylabel(r'$w^* = \sigma_{w_p}/(w_f - w_p)$', fontsize=14)

ax1.set_xlim([0,0.45])
ax1.set_xticks([0, 0.15, 0.30, 0.45])
#ax1.set_ylim([0, 0.3])
ax1.set_xticklabels([])

ax1.grid(True)

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
ax1.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
  ncol=2, borderaxespad=0)

ax1.plot(phiEval, const*phiEval**(1./3.), 'k--', zorder=1)

labelText=r"$%.1f\phi^{1/3}$" % const
ax1.text(0.15,0.15,labelText)

ax2 = fig1.add_subplot(212)
ax2.plot(phi, dwdphi2_0, 'b*', markersize=7, alpha=0.7)
ax2.plot(phi, dwdphi3_3, 'gs', markersize=7, alpha=0.7)
ax2.plot(phi, dwdphi4_0, 'ro', markersize=7, alpha=0.7)
ax2.plot(phi, dwdphi5_0, 'c^', markersize=7, alpha=0.7)

const = 0.5
ax2.plot(phiEval, 1./3.*const*phiEval**(-2./3.), 'k--', zorder=1)

ax2.set_ylim([0, 2])
ax2.set_ylabel(r"$dw^*/d\phi$")
ax2.set_xticks([0, 0.15, 0.30, 0.45])
ax2.set_xlabel(r'$\phi$', fontsize=14)

labelText=r"$\frac{1}{3} %.1f \phi^{-2/3}$" % const
ax2.text(0.15, 1.5, labelText)

ax2.grid(True)

imgname = imgdir + "hinch_stokes"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
