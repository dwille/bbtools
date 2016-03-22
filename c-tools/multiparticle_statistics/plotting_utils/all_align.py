#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
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
    alignFile = root + caseDir + '/' + 'data-tetrads/align.dat'
    nodeFile = root + caseDir + '/' + 'data-tetrads/regularNodes'

    data[stride].nnodes = file_len(nodeFile) 

    time = np.genfromtxt(alignFile, skip_header=1, usecols=0)
    data[stride].time = np.zeros(np.size(time))
    data[stride].g1s1 = np.zeros(np.size(time))
    data[stride].g2s1 = np.zeros(np.size(time))
    data[stride].g3s1 = np.zeros(np.size(time))
    data[stride].g1s2 = np.zeros(np.size(time))
    data[stride].g2s2 = np.zeros(np.size(time))
    data[stride].g3s2 = np.zeros(np.size(time))
    data[stride].g1s3 = np.zeros(np.size(time))
    data[stride].g2s3 = np.zeros(np.size(time))
    data[stride].g3s3 = np.zeros(np.size(time))
    data[stride].g1z = np.zeros(np.size(time))
    data[stride].g2z = np.zeros(np.size(time))
    data[stride].g3z = np.zeros(np.size(time))
    data[stride].s1z = np.zeros(np.size(time))
    data[stride].s2z = np.zeros(np.size(time))
    data[stride].s3z = np.zeros(np.size(time))

    data[stride].time = time - time[0]
    data[stride].g1s1 = np.genfromtxt(alignFile, skip_header=1, usecols=1)
    data[stride].g2s1 = np.genfromtxt(alignFile, skip_header=1, usecols=2)
    data[stride].g3s1 = np.genfromtxt(alignFile, skip_header=1, usecols=3)
    data[stride].g1s2 = np.genfromtxt(alignFile, skip_header=1, usecols=4)
    data[stride].g2s2 = np.genfromtxt(alignFile, skip_header=1, usecols=5)
    data[stride].g3s2 = np.genfromtxt(alignFile, skip_header=1, usecols=6)
    data[stride].g1s3 = np.genfromtxt(alignFile, skip_header=1, usecols=7)
    data[stride].g2s3 = np.genfromtxt(alignFile, skip_header=1, usecols=8)
    data[stride].g3s3 = np.genfromtxt(alignFile, skip_header=1, usecols=9)
    data[stride].g1z = np.genfromtxt(alignFile, skip_header=1, usecols=10)
    data[stride].g2z = np.genfromtxt(alignFile, skip_header=1, usecols=11)
    data[stride].g3z = np.genfromtxt(alignFile, skip_header=1, usecols=12)
    data[stride].s1z = np.genfromtxt(alignFile, skip_header=1, usecols=13)
    data[stride].s2z = np.genfromtxt(alignFile, skip_header=1, usecols=14)
    data[stride].s3z = np.genfromtxt(alignFile, skip_header=1, usecols=15)

    legendText[stride] = caseDir + ': ' + str(data[stride].nnodes)

plt.rc('font', family='serif')
color = ['r', 'g', 'b', 'k']
shade = [0.4, 0.57, 0.74, 0.9]

## Shape and Strain Alignment ##

###########
## gi_s1 ##
###########
gi_s1_Fig = plt.figure(figsize=(12,8))
gi_s1_Fig.suptitle(r'$\langle (g_i, s_1) \rangle$', fontsize=16)

# n = 500
gi_s1_500_ax = gi_s1_Fig.add_subplot(411)
for dd in range(4):
  i = dd
  gi_s1_500_ax.plot(data[i].time,data[i].g1s1, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s1_500_ax.plot(data[i].time,data[i].g2s1, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s1_500_ax.plot(data[i].time,data[i].g3s1, linewidth=3,
    color='b', alpha=shade[dd])
    
gi_s1_500_ax.set_ylabel('n = 500')
gi_s1_500_ax.set_xlim([0, 300])
gi_s1_500_ax.set_ylim([0, 0.7])
gi_s1_500_ax.grid(True)

# n = 1000
gi_s1_1000_ax = gi_s1_Fig.add_subplot(412)
for dd in range(4):
  i = 4 + dd
  gi_s1_1000_ax.plot(data[i].time,data[i].g1s1, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s1_1000_ax.plot(data[i].time,data[i].g2s1, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s1_1000_ax.plot(data[i].time,data[i].g3s1, linewidth=3,
    color='b', alpha=shade[dd])

gi_s1_1000_ax.set_ylabel('n = 1000')
gi_s1_1000_ax.set_xlim([0, 300])
gi_s1_1000_ax.set_ylim([0, 0.7])
gi_s1_1000_ax.grid(True)

# n = 1500
gi_s1_1500_ax = gi_s1_Fig.add_subplot(413)
for dd in range(4):
  i = 8 + dd
  gi_s1_1500_ax.plot(data[i].time,data[i].g1s1, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s1_1500_ax.plot(data[i].time,data[i].g2s1, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s1_1500_ax.plot(data[i].time,data[i].g3s1, linewidth=3,
    color='b', alpha=shade[dd])

gi_s1_1500_ax.set_ylabel('n = 1500')
gi_s1_1500_ax.set_xlim([0, 300])
gi_s1_1500_ax.set_ylim([0, 0.7])
gi_s1_1500_ax.grid(True)

# n = 2000
gi_s1_2000_ax = gi_s1_Fig.add_subplot(414)
for dd in range(4):
  i = 12 + dd
  gi_s1_2000_ax.plot(data[i].time,data[i].g1s1, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s1_2000_ax.plot(data[i].time,data[i].g2s1, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s1_2000_ax.plot(data[i].time,data[i].g3s1, linewidth=3,
    color='b', alpha=shade[dd])

gi_s1_2000_ax.set_ylabel('n = 2000')
gi_s1_2000_ax.set_xlim([0, 300])
gi_s1_2000_ax.set_ylim([0, 0.7])
gi_s1_2000_ax.grid(True)

plt.legend([r'$\langle(g_1, s_1)\rangle$', r'$\langle(g_2, s_1)\rangle$', r'$\langle(g_3, s_1)\rangle$'], bbox_to_anchor=(0.,-0.2,1., 0.1), 
  loc='upper left', ncol=3)

###########
## gi_s2 ##
###########
gi_s2_Fig = plt.figure(figsize=(12,8))
gi_s2_Fig.suptitle(r'$\langle (g_i, s_2) \rangle$', fontsize=16)

# n = 500
gi_s2_500_ax = gi_s2_Fig.add_subplot(411)
for dd in range(4):
  i = dd
  gi_s2_500_ax.plot(data[i].time,data[i].g1s2, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s2_500_ax.plot(data[i].time,data[i].g2s2, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s2_500_ax.plot(data[i].time,data[i].g3s2, linewidth=3,
    color='b', alpha=shade[dd])
    
gi_s2_500_ax.set_ylabel('n = 500')
gi_s2_500_ax.set_xlim([0, 300])
gi_s2_500_ax.set_ylim([0, 0.7])
gi_s2_500_ax.grid(True)

# n = 1000
gi_s2_1000_ax = gi_s2_Fig.add_subplot(412)
for dd in range(4):
  i = 4 + dd
  gi_s2_1000_ax.plot(data[i].time,data[i].g1s2, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s2_1000_ax.plot(data[i].time,data[i].g2s2, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s2_1000_ax.plot(data[i].time,data[i].g3s2, linewidth=3,
    color='b', alpha=shade[dd])

gi_s2_1000_ax.set_ylabel('n = 1000')
gi_s2_1000_ax.set_xlim([0, 300])
gi_s2_1000_ax.set_ylim([0, 0.7])
gi_s2_1000_ax.grid(True)

# n = 1500
gi_s2_1500_ax = gi_s2_Fig.add_subplot(413)
for dd in range(4):
  i = 8 + dd
  gi_s2_1500_ax.plot(data[i].time,data[i].g1s2, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s2_1500_ax.plot(data[i].time,data[i].g2s2, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s2_1500_ax.plot(data[i].time,data[i].g3s2, linewidth=3,
    color='b', alpha=shade[dd])

gi_s2_1500_ax.set_ylabel('n = 1500')
gi_s2_1500_ax.set_xlim([0, 300])
gi_s2_1500_ax.set_ylim([0, 0.7])
gi_s2_1500_ax.grid(True)

# n = 2000
gi_s2_2000_ax = gi_s2_Fig.add_subplot(414)
for dd in range(4):
  i = 12 + dd
  gi_s2_2000_ax.plot(data[i].time,data[i].g1s2, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s2_2000_ax.plot(data[i].time,data[i].g2s2, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s2_2000_ax.plot(data[i].time,data[i].g3s2, linewidth=3,
    color='b', alpha=shade[dd])

gi_s2_2000_ax.set_ylabel('n = 2000')
gi_s2_2000_ax.set_xlim([0, 300])
gi_s2_2000_ax.set_ylim([0, 0.7])
gi_s2_2000_ax.grid(True)

plt.legend([r'$\langle(g_1, s_2)\rangle$', r'$\langle(g_2, s_2)\rangle$', r'$\langle(g_3, s_2)\rangle$'], bbox_to_anchor=(0.,-0.2,1., 0.1), 
  loc='upper left', ncol=3)

###########
## gi_s3 ##
###########
gi_s3_Fig = plt.figure(figsize=(12,8))
gi_s3_Fig.suptitle(r'$\langle (g_i, s_3) \rangle$', fontsize=16)

# n = 500
gi_s3_500_ax = gi_s3_Fig.add_subplot(411)
for dd in range(4):
  i = dd
  gi_s3_500_ax.plot(data[i].time,data[i].g1s3, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s3_500_ax.plot(data[i].time,data[i].g2s3, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s3_500_ax.plot(data[i].time,data[i].g3s3, linewidth=3,
    color='b', alpha=shade[dd])
    
gi_s3_500_ax.set_ylabel('n = 500')
gi_s3_500_ax.set_xlim([0, 300])
gi_s3_500_ax.set_ylim([0, 0.7])
gi_s3_500_ax.grid(True)

# n = 1000
gi_s3_1000_ax = gi_s3_Fig.add_subplot(412)
for dd in range(4):
  i = 4 + dd
  gi_s3_1000_ax.plot(data[i].time,data[i].g1s3, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s3_1000_ax.plot(data[i].time,data[i].g2s3, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s3_1000_ax.plot(data[i].time,data[i].g3s3, linewidth=3,
    color='b', alpha=shade[dd])

gi_s3_1000_ax.set_ylabel('n = 1000')
gi_s3_1000_ax.set_xlim([0, 300])
gi_s3_1000_ax.set_ylim([0, 0.7])
gi_s3_1000_ax.grid(True)

# n = 1500
gi_s3_1500_ax = gi_s3_Fig.add_subplot(413)
for dd in range(4):
  i = 8 + dd
  gi_s3_1500_ax.plot(data[i].time,data[i].g1s3, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s3_1500_ax.plot(data[i].time,data[i].g2s3, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s3_1500_ax.plot(data[i].time,data[i].g3s3, linewidth=3,
    color='b', alpha=shade[dd])

gi_s3_1500_ax.set_ylabel('n = 1500')
gi_s3_1500_ax.set_xlim([0, 300])
gi_s3_1500_ax.set_ylim([0, 0.7])
gi_s3_1500_ax.grid(True)

# n = 2000
gi_s3_2000_ax = gi_s3_Fig.add_subplot(414)
for dd in range(4):
  i = 12 + dd
  gi_s3_2000_ax.plot(data[i].time,data[i].g1s3, linewidth=3,
    color='k', alpha=shade[dd])
  gi_s3_2000_ax.plot(data[i].time,data[i].g2s3, linewidth=3,
    color='r', alpha=shade[dd])
  gi_s3_2000_ax.plot(data[i].time,data[i].g3s3, linewidth=3,
    color='b', alpha=shade[dd])

gi_s3_2000_ax.set_ylabel('n = 2000')
gi_s3_2000_ax.set_xlim([0, 300])
gi_s3_2000_ax.set_ylim([0, 0.7])
gi_s3_2000_ax.grid(True)

plt.legend([r'$\langle(g_1, s_3)\rangle$', r'$\langle(g_2, s_3)\rangle$', r'$\langle(g_3, s_3)\rangle$'], bbox_to_anchor=(0.,-0.2,1., 0.1), 
  loc='upper left', ncol=3)

###########
## gi_z ##
###########
gi_z_Fig = plt.figure(figsize=(12,8))
gi_z_Fig.suptitle(r'$\langle (g_i, z) \rangle$', fontsize=16)

# n = 500
gi_z_500_ax = gi_z_Fig.add_subplot(411)
for dd in range(4):
  i = dd
  gi_z_500_ax.plot(data[i].time,data[i].g1z, linewidth=3,
    color='k', alpha=shade[dd])
  gi_z_500_ax.plot(data[i].time,data[i].g2z, linewidth=3,
    color='r', alpha=shade[dd])
  gi_z_500_ax.plot(data[i].time,data[i].g3z, linewidth=3,
    color='b', alpha=shade[dd])
    
gi_z_500_ax.set_ylabel('n = 500')
gi_z_500_ax.set_xlim([0, 500])
gi_z_500_ax.set_ylim([0, 0.7])
gi_z_500_ax.grid(True)

# n = 1000
gi_z_1000_ax = gi_z_Fig.add_subplot(412)
for dd in range(4):
  i = 4 + dd
  gi_z_1000_ax.plot(data[i].time,data[i].g1z, linewidth=3,
    color='k', alpha=shade[dd])
  gi_z_1000_ax.plot(data[i].time,data[i].g2z, linewidth=3,
    color='r', alpha=shade[dd])
  gi_z_1000_ax.plot(data[i].time,data[i].g3z, linewidth=3,
    color='b', alpha=shade[dd])

gi_z_1000_ax.set_ylabel('n = 1000')
gi_z_1000_ax.set_xlim([0, 500])
gi_z_1000_ax.set_ylim([0, 0.7])
gi_z_1000_ax.grid(True)

# n = 1500
gi_z_1500_ax = gi_z_Fig.add_subplot(413)
for dd in range(4):
  i = 8 + dd
  gi_z_1500_ax.plot(data[i].time,data[i].g1z, linewidth=3,
    color='k', alpha=shade[dd])
  gi_z_1500_ax.plot(data[i].time,data[i].g2z, linewidth=3,
    color='r', alpha=shade[dd])
  gi_z_1500_ax.plot(data[i].time,data[i].g3z, linewidth=3,
    color='b', alpha=shade[dd])

gi_z_1500_ax.set_ylabel('n = 1500')
gi_z_1500_ax.set_xlim([0, 500])
gi_z_1500_ax.set_ylim([0, 0.7])
gi_z_1500_ax.grid(True)

# n = 2000
gi_z_2000_ax = gi_z_Fig.add_subplot(414)
for dd in range(4):
  i = 12 + dd
  gi_z_2000_ax.plot(data[i].time,data[i].g1z, linewidth=3,
    color='k', alpha=shade[dd])
  gi_z_2000_ax.plot(data[i].time,data[i].g2z, linewidth=3,
    color='r', alpha=shade[dd])
  gi_z_2000_ax.plot(data[i].time,data[i].g3z, linewidth=3,
    color='b', alpha=shade[dd])

gi_z_2000_ax.set_ylabel('n = 2000')
gi_z_2000_ax.set_xlim([0, 500])
gi_z_2000_ax.set_ylim([0, 0.7])
gi_z_2000_ax.grid(True)

plt.legend([r'$\langle(g_1, z)\rangle$', r'$\langle(g_2, z)\rangle$', r'$\langle(g_3, s_3)\rangle$'], bbox_to_anchor=(0.,-0.2,1., 0.1), 
  loc='upper left', ncol=3)

###########
## si_z ##
###########
si_z_Fig = plt.figure(figsize=(12,8))
si_z_Fig.suptitle(r'$\langle (s_i, z) \rangle$', fontsize=16)

# n = 500
si_z_500_ax = si_z_Fig.add_subplot(411)
for dd in range(4):
  i = dd
  si_z_500_ax.plot(data[i].time,data[i].g1z, linewidth=3,
    color='k', alpha=shade[dd])
  si_z_500_ax.plot(data[i].time,data[i].g2z, linewidth=3,
    color='r', alpha=shade[dd])
  si_z_500_ax.plot(data[i].time,data[i].g3z, linewidth=3,
    color='b', alpha=shade[dd])
    
si_z_500_ax.set_ylabel('n = 500')
si_z_500_ax.set_xlim([0, 500])
si_z_500_ax.set_ylim([0, 0.7])
si_z_500_ax.grid(True)

# n = 1000
si_z_1000_ax = si_z_Fig.add_subplot(412)
for dd in range(4):
  i = 4 + dd
  si_z_1000_ax.plot(data[i].time,data[i].s1z, linewidth=3,
    color='k', alpha=shade[dd])
  si_z_1000_ax.plot(data[i].time,data[i].s2z, linewidth=3,
    color='r', alpha=shade[dd])
  si_z_1000_ax.plot(data[i].time,data[i].s3z, linewidth=3,
    color='b', alpha=shade[dd])

si_z_1000_ax.set_ylabel('n = 1000')
si_z_1000_ax.set_xlim([0, 500])
si_z_1000_ax.set_ylim([0, 0.7])
si_z_1000_ax.grid(True)

# n = 1500
si_z_1500_ax = si_z_Fig.add_subplot(413)
for dd in range(4):
  i = 8 + dd
  si_z_1500_ax.plot(data[i].time,data[i].s1z, linewidth=3,
    color='k', alpha=shade[dd])
  si_z_1500_ax.plot(data[i].time,data[i].s2z, linewidth=3,
    color='r', alpha=shade[dd])
  si_z_1500_ax.plot(data[i].time,data[i].s3z, linewidth=3,
    color='b', alpha=shade[dd])

si_z_1500_ax.set_ylabel('n = 1500')
si_z_1500_ax.set_xlim([0, 500])
si_z_1500_ax.set_ylim([0, 0.7])
si_z_1500_ax.grid(True)

# n = 2000
si_z_2000_ax = si_z_Fig.add_subplot(414)
for dd in range(4):
  i = 12 + dd
  si_z_2000_ax.plot(data[i].time,data[i].s1z, linewidth=3,
    color='k', alpha=shade[dd])
  si_z_2000_ax.plot(data[i].time,data[i].s2z, linewidth=3,
    color='r', alpha=shade[dd])
  si_z_2000_ax.plot(data[i].time,data[i].s3z, linewidth=3,
    color='b', alpha=shade[dd])

si_z_2000_ax.set_ylabel('n = 2000')
si_z_2000_ax.set_xlim([0, 500])
si_z_2000_ax.set_ylim([0, 0.7])
si_z_2000_ax.grid(True)

plt.legend([r'$\langle(s_1, z)\rangle$', r'$\langle(s_2, z)\rangle$', r'$\langle(s_3, z)\rangle$'], bbox_to_anchor=(0.,-0.2,1., 0.1), 
  loc='upper left', ncol=3)



plt.show()
