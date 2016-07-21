#!/usr/bin/env python2

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os,sys
os.system('clear')

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

print ""
print " ---- Velocity Fluctations with Volume Fraction ----"
print ""

d = 4.2
nu = .01715
rhof = 8.75e-4
partR = 2.1
g = 0.00981
rhopf = np.array([2.0, 3.3, 4.0, 5.0])
rhop = np.array([0.001750, .00289, 0.0035, 0.004375])
phi = np.array([0.087, 0.175, 0.262, 0.349])

# Setup directories
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
correlationFile = root + "vfrac_wp_correlation"
fluctFile = root + "velFlucts_raw"
termFile = root + "singlePartSedi"

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(fluctFile):
  print "      " + fluctFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(correlationFile):
  print "      " + correlationFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)
imgdir = imgdir + "f-rec-1D/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Read data
dwdphi   = np.genfromtxt(correlationFile, skip_header=1, usecols=2)
phi_sdev = np.genfromtxt(correlationFile, skip_header=1, usecols=3)
wp_sdev  = np.genfromtxt(correlationFile, skip_header=1, usecols=4)
termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)
wflucts = np.genfromtxt(fluctFile, skip_header=1, usecols=4)
wrel = np.genfromtxt(fluctFile, skip_header=1, usecols=7)

# Pull arrays of constant density, changing volume frac
dwdphi2_0 = dwdphi[0:16:4]
dwdphi3_3 = dwdphi[1:16:4]
dwdphi4_0 = dwdphi[2:16:4]
dwdphi5_0 = dwdphi[3:16:4]

phi_sdev2_0 = phi_sdev[0:16:4]
phi_sdev3_3 = phi_sdev[1:16:4]
phi_sdev4_0 = phi_sdev[2:16:4]
phi_sdev5_0 = phi_sdev[3:16:4]

wp_sdev2_0 = wp_sdev[0:16:4]
wp_sdev3_3 = wp_sdev[1:16:4]
wp_sdev4_0 = wp_sdev[2:16:4]
wp_sdev5_0 = wp_sdev[3:16:4]

# Plot dwdphi
fig1 = plt.figure(figsize=(4,4))
gs = matplotlib.gridspec.GridSpec(7,9)
ax1gs = gs[4:7,1:8]
ax2gs = gs[0:3,0:4]
ax3gs = gs[0:3,5:9]

ax1 = fig1.add_subplot(ax1gs)
ax1.plot(phi, dwdphi2_0, 'ks', color='0.0')
ax1.plot(phi, dwdphi3_3, 'ko', color='0.2')
ax1.plot(phi[0:3], dwdphi4_0[0:3], 'k^', color='0.4')
ax1.plot(phi[0:2], dwdphi5_0[0:2], 'kx', color='0.6')

ax1.set_xlabel(r"$\phi$")
ax1.set_xlim([0, 0.45])
ax1.set_ylabel(r"$\frac{dw_p^\prime}{d\phi}$")

ax1.xaxis.set_major_locator(MultipleLocator(0.15))
ax1.xaxis.set_minor_locator(MultipleLocator(0.075))
ax1.yaxis.set_major_locator(MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))

# Plot sdev of phi

ax2 = fig1.add_subplot(ax2gs)
ax2.plot(phi, phi_sdev2_0, 'ks', color='0.0')
ax2.plot(phi, phi_sdev3_3, 'ko', color='0.2')
ax2.plot(phi[0:3], phi_sdev4_0[0:3], 'k^', color='0.4')
ax2.plot(phi[0:2], phi_sdev5_0[0:2], 'kx', color='0.6')

ax2.set_xlabel(r"$\phi$")
ax2.set_xlim([0, 0.45])
ax2.set_ylabel(r"$\phi_{sdev}^\prime$")
ax2.set_ylim([0, 0.03])

ax2.xaxis.set_major_locator(MultipleLocator(0.15))
ax2.xaxis.set_minor_locator(MultipleLocator(0.075))
ax2.yaxis.set_major_locator(MultipleLocator(0.01))
ax2.yaxis.set_minor_locator(MultipleLocator(0.005))

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
ax2.legend(lText, bbox_to_anchor=(.25,1.10,1.8,1),
  loc="lower left", mode="expand", ncol=2, borderaxespad=0)

# Plot sdev of wp

ax3 = fig1.add_subplot(ax3gs)
ax3.plot(phi, wp_sdev2_0, 'ks', color='0.0')
ax3.plot(phi, wp_sdev3_3, 'ko', color='0.2')
ax3.plot(phi[0:3], wp_sdev4_0[0:3], 'k^', color='0.4')
ax3.plot(phi[0:2], wp_sdev5_0[0:2], 'kx', color='0.6')

ax3.set_xlabel(r"$\phi$")
ax3.set_xlim([0, 0.45])
ax3.set_ylabel(r"$w_{p,sdev}^\prime$")
ax3.set_ylim([0, 0.025])

ax3.xaxis.set_major_locator(MultipleLocator(0.15))
ax3.xaxis.set_minor_locator(MultipleLocator(0.075))
ax3.yaxis.set_major_locator(MultipleLocator(0.01))
ax3.yaxis.set_minor_locator(MultipleLocator(0.005))

ax3.yaxis.tick_right()
ax3.yaxis.set_ticks_position('both')
ax3.yaxis.set_label_position('right')

imgname = imgdir + "vfrac_wp_correlations"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
