#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os, sys
os.system('clear')

print ""
print " ---- Wavespeed Plotting Utility ---- "
print ""

d = 4.2
nu = .01715

root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/simdata/"
waveFile = root + "waveData"
termFile = root + "singlePartSedi"

print "      Root directory set to: " + root
print "      Using location: " + waveFile
print ""
print "      Using d = %.2f" % d
print "      Using nu = %.2f" % nu

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(waveFile):
  print "      " + waveFile + " does not exist. Exiting..."
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
nVal = np.genfromtxt(waveFile, skip_header=1, usecols=0)
rhoVal = np.genfromtxt(waveFile, skip_header=1, usecols=1)
waveSpeed = np.genfromtxt(waveFile, skip_header=1, usecols=2)
termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)

n = np.unique(nVal)
phi = np.array([0.087, 0.175, 0.262, 0.349])
rho = np.unique(rhoVal)

# arrays-- constant density, variable vfrac
rho20 = np.array([waveSpeed[0], waveSpeed[4], waveSpeed[8], waveSpeed[12]])
rho33 = np.array([waveSpeed[1], waveSpeed[5], waveSpeed[9], waveSpeed[13]])
rho40 = np.array([waveSpeed[2], waveSpeed[6], waveSpeed[10]])
rho50 = np.array([waveSpeed[3], waveSpeed[7]])

# Garside and Al-Dibouni n coeff
tmp = 0.1*Ret**0.9
nGAD = (5.1 + 2.7*tmp)/(1 + tmp)

# Wallis relations
Ret = d*termVel/nu

phiEval = np.linspace(0.05,0.6,101)
wallis_speed = np.zeros((np.size(phiEval),4))
for rr,_ in enumerate(rho):
  for pp,_ in enumerate(phiEval):
    wallis_speed[pp,rr] = nGAD[rr]*(1-phiEval[pp])**(nGAD[rr] - 1)*phiEval[pp]*termVel[rr] * 1000


##
 # Plots
 ##

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

# Data
ax1.plot(phi, rho20, 'b*', alpha=0.9,markersize=6)
ax1.plot(phi, rho33, 'gs', alpha=0.9,markersize=6)
ax1.plot(phi[0:3], rho40, 'ro', alpha=0.9,markersize=6)
ax1.plot(phi[0:2], rho50, 'c^', alpha=0.9,markersize=6)

ax1.set_xlabel(r'$\phi$', fontsize=14)
ax1.set_ylabel(r'$c\ [mm/s]$', fontsize=14)
ax1.set_xlim([0,0.5])
ax1.set_xticks([0, 0.10, 0.20, 0.30, 0.40, 0.50])
ax1.grid(True)

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
ax1.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
  ncol=2, borderaxespad=0)

# Wallis correlations
ax1.plot(phiEval, wallis_speed[:,0], 'b--',zorder=1)
ax1.plot(phiEval, wallis_speed[:,1], 'g--',zorder=1)
ax1.plot(phiEval, wallis_speed[:,2], 'r--',zorder=1)
ax1.plot(phiEval, wallis_speed[:,3], 'c--',zorder=1)

imgname = imgdir + "wavespeed"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      ... Done!"
