#!/usr/bin/env python2

import sys, os
import glob
import matplotlib.pyplot as plt
from matplotlib import lines as mlines
import numpy as np
import re

os.system('clear')

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line

print ""
print " ---- Anisotropy Measures Plotting Utility ---- "
print "                Kurtosis"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "      Sim root directory set to: " + root

# Parameter sweep
partList = ['500', '1000', '1500', '2000']
densList = ['rho2.0', 'rho3.3', 'rho4.0', 'rho5.0']
legendText = ['']*16

# Initialize arrays
nTetrads = np.zeros(16)
data = [ structtype() for i in range(16) ]
RoG = [ structtype() for i in range(16) ]
EVar = [ structtype() for i in range(16) ]
Shape = [ structtype() for i in range(16) ]
I1 = [ structtype() for i in range(16) ]
I2 = [ structtype() for i in range(16) ]
I3 = [ structtype() for i in range(16) ]
 
# Loop over all directory's
for pp, part in enumerate(partList):
  for dd, dens in enumerate(densList):
    stride = 4*pp + dd
    caseDir = part + '/' + dens
    infoFile = root + caseDir + '/data-tetrads/info.dat'
    nodeFile = root + caseDir + '/data-tetrads/regularNodes'
    statkurt = root + caseDir + '/data-tetrads/stat.kurt'

    nTetrads[stride] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

    tmpT = np.genfromtxt(statkurt, skip_header=1, usecols=0)
    data[stride].time = np.zeros(np.size(tmpT))
    data[stride].time = tmpT - tmpT[0]

    # kurt
    RoG[stride].kurt   = np.zeros(np.size(tmpT))
    EVar[stride].kurt  = np.zeros(np.size(tmpT))
    Shape[stride].kurt = np.zeros(np.size(tmpT))
    I1[stride].kurt    = np.zeros(np.size(tmpT))
    I2[stride].kurt    = np.zeros(np.size(tmpT))
    I3[stride].kurt    = np.zeros(np.size(tmpT))

    RoG[stride].kurt   = np.genfromtxt(statkurt, skip_header=1, usecols=1)
    EVar[stride].kurt  = np.genfromtxt(statkurt, skip_header=1, usecols=2)
    Shape[stride].kurt = np.genfromtxt(statkurt, skip_header=1, usecols=3)
    I1[stride].kurt    = np.genfromtxt(statkurt, skip_header=1, usecols=4)
    I2[stride].kurt    = np.genfromtxt(statkurt, skip_header=1, usecols=5)
    I3[stride].kurt    = np.genfromtxt(statkurt, skip_header=1, usecols=6)

    legendText[stride] = caseDir + ': ' + str(nTetrads[stride])

# Plotting specs
plt.rc('font', family='serif')
colors = ['r', 'g', 'b', 'k']
shades = [0.4, 0.57, 0.74, 0.9]
fsAx = 14
fSize = (12,8)
lWidth = 2

# Legend specs
rho2_spec = mlines.Line2D([],[], color='k', alpha=shades[0], linewidth=lWidth,
 label=r'$\rho^* = 2.0$')
rho3_spec = mlines.Line2D([],[], color='k', alpha=shades[1], linewidth=lWidth,
 label=r'$\rho^* = 3.3$')
rho4_spec = mlines.Line2D([],[], color='k', alpha=shades[2], linewidth=lWidth,
 label=r'$\rho^* = 4.0$')
rho5_spec = mlines.Line2D([],[], color='k', alpha=shades[3], linewidth=lWidth,
 label=r'$\rho^* = 5.0$')

g1si_spec = mlines.Line2D([],[], color='k', linewidth=lWidth,
  label=r'$\langle(g_1, s_i)\rangle$')
g2si_spec = mlines.Line2D([],[], color='r', linewidth=lWidth,
  label=r'$\langle(g_2, s_i)\rangle$')
g3si_spec = mlines.Line2D([],[], color='b', linewidth=lWidth,
  label=r'$\langle(g_3, s_i)\rangle$')

## Radius of Gyration
rgFig = plt.figure(figsize=(12,8))
rgFig.suptitle('Kurtosis of Radius of Gyration', fontsize=16)
rg_ax = rgFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    rg_ax.semilogx(data[i].time, RoG[i].kurt, linewidth=3,
      color=colors[pp], alpha=shades[dd])

#xpnts = np.array([100, 10000])
#ypnts = np.power(xpnts, 0.6) / 5.
#rg_ax.loglog(xpnts, ypnts, 'k--', linewidth=3)
#rg_ax.text(300, 5, 'slope = 0.6')

rg_ax.legend(legendText, ncol=2,loc='upper left',framealpha=0.8)
rg_ax.set_xlabel("Time [ms]", fontsize=16)
rg_ax.set_ylabel(r"$R$ kurtosis", fontsize=16)

# EVar
vFig = plt.figure(figsize=(12,8))
vFig.suptitle('Kurtosis of Eigenvalue variance', fontsize=16)
v_ax = vFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    v_ax.semilogx(data[i].time, EVar[i].kurt, linewidth=3,
      color=colors[pp], alpha=shades[dd])

v_ax.set_xlabel("Time [ms]", fontsize=16)
v_ax.set_ylabel("Eigenvalue Variance - kurtosis", fontsize=16)
v_ax.tick_params(which='major', length=10)
v_ax.tick_params(which='minor', length=7)

#xpnts = np.array([10, 90])
#ypnts = np.power(xpnts, 0.7) / 37
#v_ax.loglog(xpnts, ypnts, 'k--', linewidth=3)
#v_ax.text(10, 0.4, 'slope = 0.7')


# Shape
sFig = plt.figure(figsize=(12,8))
sFig.suptitle('Kurtosis of Shape', fontsize=16)
s_ax = sFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    s_ax.semilogx(data[i].time, Shape[i].kurt, linewidth=3,
      color=colors[pp], alpha=shades[dd])

s_ax.set_xlabel("Time [ms]", fontsize=16)
s_ax.set_ylabel("Shape kurtosis", fontsize=16)
s_ax.legend(legendText, ncol=2,loc='upper right',framealpha=0.8)

#xpnts = np.array([9, 110])
#ypnts = np.power(xpnts, 1.75) / 2000
#s_ax.loglog(xpnts, ypnts, 'k--', linewidth=3)
#s_ax.text(10, 0.4, 'slope = 1.75')

## I_j ##
iFig = plt.figure(figsize=(12,8))
iFig.suptitle('Kurtosis of I', fontsize=16)
i_ax = iFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    i_ax.semilogx(data[i].time, I1[i].kurt, linewidth=3,
      color=colors[pp], alpha=shades[dd])
    i_ax.semilogx(data[i].time, I2[i].kurt, linewidth=3,
      color=colors[pp], alpha=shades[dd])
    i_ax.semilogx(data[i].time, I3[i].kurt, linewidth=3,
      color=colors[pp], alpha=shades[dd])

i_ax.text(5, 0.55, 'I1', fontsize=16)
i_ax.text(5, 0.375, 'I2', fontsize=16)
i_ax.text(5, 0.175, 'I3', fontsize=16)
i_ax.set_ylabel("I kurtosis", fontsize=16)

plt.tight_layout()
plt.show()
