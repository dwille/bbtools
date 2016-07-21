#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import matplotlib.lines as mlines
os.system('clear')

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

print ""
print " ---- Wavespeed Plotting Utility ---- "
print ""

d = 4.2
nu = .01715
g = 0.00981
rho_f = 8.75e-4
rho_p = np.array([0.001750, 0.00289, 0.0035, 0.004375])

home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
phaseFile = root + "phaseAveragedFluidVel"
waveFile = root + "waveData"
termFile = root + "singlePartSedi"
rzFile = root + "richardson_zaki_coeffs"

print "      Root directory set to: " + root
print "      Using location: " + waveFile
print ""
print "      Using d = %.2f" % d
print "      Using nu = %.4f" % nu

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(waveFile):
  print "      " + waveFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(phaseFile):
  print "      " + phaseFile + " does not exist. Exiting..."
  print ""
  sys.exit()
 
# Create imgdir if necessary
imgdir = root + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

imgdir = imgdir + "f-rec-1D/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Read Data
npartsVal = np.genfromtxt(waveFile, skip_header=1, usecols=0)
rhoVal = np.genfromtxt(waveFile, skip_header=1, usecols=1)
# wavespeed in mm/s
waveSpeed = np.genfromtxt(waveFile, skip_header=1, usecols=2)
# termVel and flowVel in mm/ms
termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)
flowVel = np.genfromtxt(phaseFile, skip_header=1, usecols=1)
n = np.genfromtxt(rzFile, skip_header=1, usecols=1)
kappa = np.genfromtxt(rzFile, skip_header=1, usecols=2)

# convert all to m/s
# in mm/s, divide by 1000
# in mm/ms, convert to m/s by *1
waveSpeed /= 1000
termVel *= 1.
flowVel *= 1.

nParts = np.unique(npartsVal)
phi = np.array([0.087, 0.175, 0.262, 0.349])
rho = np.unique(rhoVal)

# wavespeed -- constant density, variable vfrac
rho20_c = waveSpeed[0:16:4] / termVel[0]
rho33_c = waveSpeed[1:16:4] / termVel[1]
rho40_c = waveSpeed[2:12:4] / termVel[2]
rho50_c = waveSpeed[3:8:4]  / termVel[3]

# fluid velocity -- constant density, variable vfrac
rho20_vel = flowVel[0:16:4]
rho33_vel = flowVel[1:16:4]
rho40_vel = flowVel[2:12:4]
rho50_vel = flowVel[3:8:4]

# Wallis relations
phiEval = np.linspace(0.05,0.6,101)
wallis_speed = np.zeros((np.size(phiEval),4))
#gibi_speed = np.zeros((np.size(phiEval),4))
for rr,_ in enumerate(rho):
  for pp,_ in enumerate(phiEval):
    wallis_speed[pp,rr] = termVel[rr]*kappa[rr]*n[rr]*phiEval[pp]*(1-phiEval[pp])**(n[rr] - 1)
    wallis_speed[pp,rr] = kappa[rr]*n[rr]*phiEval[pp]*(1-phiEval[pp])**(n[rr] - 1)

##
 # Wave Speed Plots
 ##
fig1 = plt.figure(figsize=(5,2.5))
ax1 = fig1.add_subplot(111)

# dark
bl= "#4C72B0"
gr= "#55A868"
re= "#C44E52"
pu= "#8172B2"
go= "#CCB974"
cy= "#64B5CD"

# Data
ax1.plot(phi, rho20_c, '*', color=bl, markersize=7)
ax1.plot(phi, rho33_c, 's', color=gr, markersize=7)
ax1.plot(phi[0:3], rho40_c, 'o', color=re, markersize=7)
ax1.plot(phi[0:2], rho50_c, '^', color=cy, markersize=7)

ax1.set_xlabel(r'$\phi$')
#ax1.set_ylabel(r'$c$')
ax1.set_ylabel(r'$c/w_t$')
ax1.set_xlim([0,0.5])
ax1.grid(True)

ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
h1 = mlines.Line2D([],[], linestyle=':', color=bl, marker='*', label=lText[0])
h2 = mlines.Line2D([],[], linestyle='-.', color=gr, marker='s', label=lText[1])
h3 = mlines.Line2D([],[], linestyle='-', color=re, marker='o', label=lText[2])
h4 = mlines.Line2D([],[], linestyle='--', color=cy, marker='^', label=lText[3])
ax1.legend(handles=[h1,h2,h3,h4], bbox_to_anchor=(0,1.05,1,1), loc="lower left",
  mode="expand", ncol=2, borderaxespad=0)

# Wallis correlations
ax1.plot(phiEval, wallis_speed[:,0], ':', color=bl, zorder=1)
ax1.plot(phiEval, wallis_speed[:,1], '-.', color=gr, zorder=1)
ax1.plot(phiEval, wallis_speed[:,2], '-', color=re, zorder=1)
ax1.plot(phiEval, wallis_speed[:,3], '--', color=cy, zorder=1)


#imgname = imgdir + "wavespeed"
imgname = imgdir + "wavespeed_normed"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# ##
#  # Froude Number Plots
#  ##
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
## Froude number
#rho20_Fr = rho20_vel/rho20_c
#rho33_Fr = rho33_vel/rho33_c
#rho40_Fr = rho40_vel/rho40_c
#rho50_Fr = rho50_vel/rho50_c
#
# 
# ax2.plot(phi, rho20_Fr, 'b*', alpha=0.9,markersize=6)
# ax2.plot(phi, rho33_Fr, 'gs', alpha=0.9,markersize=6)
# ax2.plot(phi[0:3], rho40_Fr, 'ro', alpha=0.9,markersize=6)
# ax2.plot(phi[0:2], rho50_Fr, 'c^', alpha=0.9,markersize=6)
# 
# ax2.set_xlabel(r'$\phi$', fontsize=14)
# ax2.set_ylabel(r'$Fr = \frac{w_f}{c}$', fontsize=14)
# ax2.set_xlim([0,0.5])
# ax2.set_xticks([0, 0.10, 0.20, 0.30, 0.40, 0.50])
# ax2.grid(True)
# 
# lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
#          r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
# ax2.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
#   ncol=2, borderaxespad=0)
# 
# imgname = imgdir + "froude"
# plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
# plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      ... Done!"
