#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import re

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line

print "Anisotropy Measures Plotting Utility"
print "   Plot measures for all simulations"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print " Sim root directory set to: " + root

# Parameter sweep
partList = ['500', '1000', '1500', '2000']
densList = ['rho2.0', 'rho3.3', 'rho4.0', 'rho5.0']
legendText = ['']*16

# Initialize mean/std arrays
data = [ structtype() for i in range(16) ]
 
# Loop over all directory's
for pp, part in enumerate(partList):
  for dd, dens in enumerate(densList):
    stride = 4*pp + dd
    caseDir = part + '/' + dens
    statFile = root + caseDir + '/' + 'data-tetrads/stat.dat'
    nodeFile = root + caseDir + '/' + 'data-tetrads/regularNodes'

    data[stride].nnodes = file_len(nodeFile) 

    time = np.genfromtxt(statFile, skip_header=1, usecols=0)
    data[stride].time = np.zeros(np.size(time))
    data[stride].meanR2 = np.zeros(np.size(time))
    data[stride].meanVar = np.zeros(np.size(time))
    data[stride].meanShape = np.zeros(np.size(time))

    data[stride].time = time - time[0]
    data[stride].meanR2 = np.genfromtxt(statFile, skip_header=1, usecols=1)
    data[stride].meanVar = np.genfromtxt(statFile, skip_header=1, usecols=2)
    data[stride].meanShape = np.genfromtxt(statFile, skip_header=1, usecols=3)

    legendText[stride] = caseDir + ': ' + str(data[stride].nnodes)

plt.rc('font', family='serif')
colors = ['r', 'g', 'b', 'k']
shades = [0.4, 0.57, 0.74, 0.9]

# Radius of Gyration
rgFig = plt.figure(figsize=(12,8))
rgFig.suptitle('Radius of Gyration', fontsize=16)
rg_ax = rgFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    rg_ax.loglog(data[i].time, np.sqrt(data[i].meanR2), linewidth=3,
      color=colors[pp], alpha=shades[dd])

xpnts = np.array([100, 10000])
ypnts = np.power(xpnts, 1) / 10
rg_ax.loglog(xpnts, ypnts, 'k--', linewidth=3)
rg_ax.text(1000, 20, 'slope = 1')
rg_ax.legend(legendText, ncol=2,loc='upper left')
rg_ax.set_xlabel("Time [ms]", fontsize=16)
rg_ax.set_ylabel(r"$\langle R^2 \rangle$", fontsize=16)

# Var
vFig = plt.figure(figsize=(12,8))
vFig.suptitle('Shape', fontsize=16)
v_ax = vFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    v_ax.semilogx(data[i].time, data[i].meanVar, linewidth=3,
      color=colors[pp], alpha=shades[dd])

v_ax.set_xlabel("Time [ms]", fontsize=16)
v_ax.set_ylabel("Variance", fontsize=16)

# Shape
sFig = plt.figure(figsize=(12,8))
sFig.suptitle('Shape', fontsize=16)
s_ax = sFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    s_ax.semilogx(data[i].time, data[i].meanShape, linewidth=3,
      color=colors[pp], alpha=shades[dd])

s_ax.set_xlabel("Time [ms]", fontsize=16)
s_ax.set_ylabel("Shape", fontsize=16)

plt.tight_layout()
plt.show()
