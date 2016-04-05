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
print "                   Skewness"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "      Sim root directory set to: " + root

# Create imgdir if necessary
imgdir = root + "img/shape/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

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
    statskew = root + caseDir + '/data-tetrads/stat.skew'

    nTetrads[stride] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

    tmpT = np.genfromtxt(statskew, skip_header=1, usecols=0)
    data[stride].time = np.zeros(np.size(tmpT))
    data[stride].time = tmpT - tmpT[0]

    # skew
    RoG[stride].skew   = np.zeros(np.size(tmpT))
    EVar[stride].skew  = np.zeros(np.size(tmpT))
    Shape[stride].skew = np.zeros(np.size(tmpT))
    I1[stride].skew    = np.zeros(np.size(tmpT))
    I2[stride].skew    = np.zeros(np.size(tmpT))
    I3[stride].skew    = np.zeros(np.size(tmpT))

    RoG[stride].skew   = np.genfromtxt(statskew, skip_header=1, usecols=1)
    EVar[stride].skew  = np.genfromtxt(statskew, skip_header=1, usecols=2)
    Shape[stride].skew = np.genfromtxt(statskew, skip_header=1, usecols=3)
    I1[stride].skew    = np.genfromtxt(statskew, skip_header=1, usecols=4)
    I2[stride].skew    = np.genfromtxt(statskew, skip_header=1, usecols=5)
    I3[stride].skew    = np.genfromtxt(statskew, skip_header=1, usecols=6)

    legendText[stride] = caseDir + ': ' + str(nTetrads[stride])

# Plot specs
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=10, numpoints=3)
plt.rc('lines', markersize=4, linewidth=2)
labelx = -0.17
colors = ['r', 'g', 'b', 'k']
shades = [0.4, 0.57, 0.74, 0.9]

## Radius of Gyration
rgFig = plt.figure()
rg_ax = rgFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    rg_ax.semilogx(data[i].time,RoG[i].skew, color=colors[pp], alpha=shades[dd])

rg_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))
rg_ax.set_xlabel("Time [ms]")
rg_ax.set_ylabel(r"$\mathrm{Skew}[R_g]$")
# Save
imgname = imgdir + "all_rog_skew"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# EVar
vFig = plt.figure()
v_ax = vFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    v_ax.semilogx(data[i].time, EVar[i].skew, color=colors[pp], 
      alpha=shades[dd])

v_ax.set_xlabel("Time [ms]")
v_ax.set_ylabel(r"$\mathrm{Skew}[\Delta]$")
v_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# Save
imgname = imgdir + "all_evar_skew"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Shape
sFig = plt.figure()
s_ax = sFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    s_ax.semilogx(data[i].time, Shape[i].skew, color=colors[pp], 
      alpha=shades[dd])

s_ax.set_xlabel("Time [ms]")
s_ax.set_ylabel(r"$\mathrm{Skew}[S]$")
s_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# Save
imgname = imgdir + "all_shape_skew"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## I_j ##
iFig = plt.figure()
i_ax = iFig.add_subplot(111)
for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    i_ax.semilogx(data[i].time, I1[i].skew, color=colors[pp],
      alpha=shades[dd])
    i_ax.semilogx(data[i].time, I2[i].skew, color=colors[pp],
      alpha=shades[dd], label='_nolegend_')
    i_ax.semilogx(data[i].time, I3[i].skew, color=colors[pp],
      alpha=shades[dd], label='_nolegend_')

i_ax.text(1.25, -0.50, 'I1')
i_ax.text(1.25, 0.00, 'I2')
i_ax.text(1.25, 0.50, 'I3')
i_ax.set_ylabel(r"$\mathrm{Skew}[I_j]$")

i_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# Save
imgname = imgdir + "all_i_skew"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
