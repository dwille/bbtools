#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os
os.system('clear')

## GET INFO
print "   Phase-Averaged Fluid Velocity Plotting Utility"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "Root directory set to: " + root
dataFile = root + "phaseAveragedFluidVel"
print ""
print "Using location: " + dataFile

# Read data
flowVel = np.genfromtxt(dataFile, skip_header=1, usecols=1)

# Constant rho -- alpha increasing
rho2_0 = flowVel[0:16:4]
rho3_3 = flowVel[1:16:4]
rho4_0 = flowVel[2:16:4]
rho5_0 = flowVel[3:16:4]

# Constant alpha -- rho increasing
n0500 = flowVel[0:4]
n1000 = flowVel[4:8]
n1500 = flowVel[8:12]
n2000 = flowVel[12:16]

rho = np.array([2.0, 3.3, 4.0, 5.0])
alpha = np.array([0.087, 0.175, 0.262, 0.349])

# Plot
flow = plt.figure(figsize=(12,8))
flow.suptitle('Phase Averaged Fluid Velocity', fontsize=20)

# Constant rho -- alpha increasing
nAx = flow.add_subplot(211)
nAx.plot(alpha,rho2_0, 'o--', linewidth=3, markersize=10)
nAx.plot(alpha,rho3_3, 'o--', linewidth=3, markersize=10)
nAx.plot(alpha,rho4_0, 'o--', linewidth=3, markersize=10)
nAx.plot(alpha,rho5_0, 'o--', linewidth=3, markersize=10)

nAx.set_xlabel('Particle Volume Fraction', fontsize=16)
nAx.set_ylabel('W_f', fontsize=16)

nAx.set_xlim([0,0.4])
nAx.set_ylim([0,0.4])

nAx.grid(True)

nAx.legend(['rho* = 2.0', 'rho* = 3.3', 'rho* = 4.0', 'rho* = 5.0'])

# Constant alpha -- rho increasing
rhoAx = flow.add_subplot(212)
rhoAx.plot(rho, n0500, 'o--', linewidth=3, markersize=10)
rhoAx.plot(rho, n1000, 'o--', linewidth=3, markersize=10)
rhoAx.plot(rho, n1500, 'o--', linewidth=3, markersize=10)
rhoAx.plot(rho, n2000, 'o--', linewidth=3, markersize=10)

rhoAx.set_xlabel('Density Ratio', fontsize=16)
rhoAx.set_ylabel('W_f', fontsize=16)

rhoAx.set_xlim([1,6])
rhoAx.set_ylim([0,0.4])

rhoAx.grid(True)

rhoAx.legend(['alpha = 8.7%', 'alpha = 17.5%', 'alpha = 26.2%', 'alpha = 34.9%'],
  loc='lower right')



plt.show()
