#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import re

## Get info
print ""
print " ---- Alignment plotting utility ---- "
print ""
root = '/home-1/dwillen3@jhu.edu/scratch/triply_per/'
simdir = raw_input("      Simulation root: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'
datadir = root + simdir + 'data-tetrads/'
print "      Data directory set to: " + datadir
alignFile = datadir + 'align.mean'

time = np.genfromtxt(alignFile, skip_header=1, usecols=0)
time = time - time[0]

g1s1 = np.genfromtxt(alignFile, skip_header=1, usecols=1)
g2s1 = np.genfromtxt(alignFile, skip_header=1, usecols=2)
g3s1 = np.genfromtxt(alignFile, skip_header=1, usecols=3)

g1s2 = np.genfromtxt(alignFile, skip_header=1, usecols=4)
g2s2 = np.genfromtxt(alignFile, skip_header=1, usecols=5)
g3s2 = np.genfromtxt(alignFile, skip_header=1, usecols=6)

g1s3 = np.genfromtxt(alignFile, skip_header=1, usecols=7)
g2s3 = np.genfromtxt(alignFile, skip_header=1, usecols=8)
g3s3 = np.genfromtxt(alignFile, skip_header=1, usecols=9)

g1z = np.genfromtxt(alignFile, skip_header=1, usecols=10)
g2z = np.genfromtxt(alignFile, skip_header=1, usecols=11)
g3z = np.genfromtxt(alignFile, skip_header=1, usecols=12)

s1z = np.genfromtxt(alignFile, skip_header=1, usecols=13)
s2z = np.genfromtxt(alignFile, skip_header=1, usecols=14)
s3z = np.genfromtxt(alignFile, skip_header=1, usecols=15)

wz = np.genfromtxt(alignFile, skip_header=1, usecols=16)

wg1 = np.genfromtxt(alignFile, skip_header=1, usecols=17)
wg2 = np.genfromtxt(alignFile, skip_header=1, usecols=18)
wg3 = np.genfromtxt(alignFile, skip_header=1, usecols=19)

ws1 = np.genfromtxt(alignFile, skip_header=1, usecols=20)
ws2 = np.genfromtxt(alignFile, skip_header=1, usecols=21)
ws3 = np.genfromtxt(alignFile, skip_header=1, usecols=22)

#wnorm = np.genfromtxt(alignFile, skip_header=1, usecols=23)

## Principal axes / strain ##
plt.rc('font', family='serif')
g_s = plt.figure(figsize=(12,8))

g_s.suptitle('Alignment of Shape Principal Axes with Strain', fontsize=16)

# major
g1s_Ax = g_s.add_subplot(311)
g1s_Ax.plot(time, g1s1, 'ko-', linewidth=1.5, markevery=1)
g1s_Ax.plot(time, g1s2, 'bo-', linewidth=1.5, markevery=1)
g1s_Ax.plot(time, g1s3, 'ro-', linewidth=1.5, markevery=1)

g1s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_1 \cdot \mathbf{s(0)}_i)$', 
  fontsize=15)
#g1s_Ax.set_ylim([0, 0.6])
g1s_Ax.set_xlim([0, 200])
g1s_Ax.xaxis.set_minor_locator(AutoMinorLocator())
g1s_Ax.yaxis.set_minor_locator(AutoMinorLocator())
g1s_Ax.tick_params(which='major', length=6)
g1s_Ax.tick_params(which='minor', length=3)

# middle
g2s_Ax = g_s.add_subplot(312)
g2s_Ax.plot(time, g2s1, 'ko-', linewidth=1.5, markevery=1)
g2s_Ax.plot(time, g2s2, 'bo-', linewidth=1.5, markevery=1)
g2s_Ax.plot(time, g2s3, 'ro-', linewidth=1.5, markevery=1)


g2s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_2 \cdot \mathbf{s(0)}_i)$', 
  fontsize=15)
#g2s_Ax.set_ylim([0, 0.6])
g2s_Ax.set_xlim([0, 200])
g2s_Ax.xaxis.set_minor_locator(AutoMinorLocator())
g2s_Ax.yaxis.set_minor_locator(AutoMinorLocator())
g2s_Ax.tick_params(which='major', length=6)
g2s_Ax.tick_params(which='minor', length=3)

g3s_Ax = g_s.add_subplot(313)
g3s_Ax.plot(time, g3s1, 'ko-', linewidth=1.5, markevery=10)
g3s_Ax.plot(time, g3s2, 'bo-', linewidth=1.5, markevery=10)
g3s_Ax.plot(time, g3s3, 'ro-', linewidth=1.5, markevery=10)

# minor
g3s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_3 \cdot \mathbf{s(0)}_i)$', 
  fontsize=15)
#g3s_Ax.set_ylim([0, 0.6])
g3s_Ax.set_xlim([0, 200])
g3s_Ax.set_xlabel('Time [ms]')
g3s_Ax.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])
g3s_Ax.xaxis.set_minor_locator(AutoMinorLocator())
g3s_Ax.yaxis.set_minor_locator(AutoMinorLocator())
g3s_Ax.tick_params(which='major', length=6)
g3s_Ax.tick_params(which='minor', length=3)

plt.show()
## vorticity / paxes, strain ##
wFig = plt.figure(figsize=(12,8))
wFig.suptitle('Alignment of shape and strain with vorticity', fontsize=16)

# principal axes
gw_Ax = wFig.add_subplot(311)
gw_Ax.plot(time, wg1, 'ko-', linewidth=1.5, markevery=10)
gw_Ax.plot(time, wg2, 'bo-', linewidth=1.5, markevery=10)
gw_Ax.plot(time, wg3, 'ro-', linewidth=1.5, markevery=10)

gw_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g(0)}_i \cdot \mathbf{\omega})$', 
  fontsize=15)
gw_Ax.legend(['Major Shape Axis', 'Middle Shape Axis', 'Minor Shape Axis'])
gw_Ax.xaxis.set_minor_locator(AutoMinorLocator())
gw_Ax.yaxis.set_minor_locator(AutoMinorLocator())
gw_Ax.tick_params(which='major', length=6)
gw_Ax.tick_params(which='minor', length=3)

# strain
sw_Ax = wFig.add_subplot(312)
sw_Ax.plot(time, ws1, 'ko-', linewidth=1.5, markevery=10)
sw_Ax.plot(time, ws2, 'bo-', linewidth=1.5, markevery=10)
sw_Ax.plot(time, ws3, 'ro-', linewidth=1.5, markevery=10)
sw_Ax.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])

sw_Ax.set_ylabel(r'$\cos \theta = (\mathbf{s(0)}_i \cdot \mathbf{\omega})$', 
  fontsize=15)
sw_Ax.xaxis.set_minor_locator(AutoMinorLocator())
sw_Ax.yaxis.set_minor_locator(AutoMinorLocator())
sw_Ax.tick_params(which='major', length=6)
sw_Ax.tick_params(which='minor', length=3)

# vorticity magnitude
w_Ax = wFig.add_subplot(313)

#w_Ax.plot(time, wnorm, 'ko-', linewidth=1.5)
#
#w_Ax.set_xlabel('Time [ms]')
#w_Ax.set_ylabel('Vorticity')
#w_Ax.xaxis.set_minor_locator(AutoMinorLocator())
#w_Ax.yaxis.set_minor_locator(AutoMinorLocator())
#w_Ax.tick_params(which='major', length=6)
#w_Ax.tick_params(which='minor', length=3)

## alignment with z ##
g_z = plt.figure(figsize=(12,8))
g_z.suptitle('Alignment of shape, strain, and vorticity with gravity', 
  fontsize=16)

gzAx = g_z.add_subplot(311)
gzAx.plot(time, g1z, 'ko-', linewidth=1.5, markevery=10)
gzAx.plot(time, g2z, 'bo-', linewidth=1.5, markevery=10)
gzAx.plot(time, g3z, 'ro-', linewidth=1.5, markevery=10)

gzAx.set_ylim([0, 1])
gzAx.set_ylabel(r'$\cos \theta = (\mathbf{g}_i \cdot \mathbf{z})$', fontsize=15)
gzAx.legend(['Major Shape Axis', 'Middle Shape Axis', 'Minor Shape Axis'],
  loc='upper left')
gzAx.xaxis.set_minor_locator(AutoMinorLocator())
gzAx.yaxis.set_minor_locator(AutoMinorLocator())
gzAx.tick_params(which='major', length=6)
gzAx.tick_params(which='minor', length=3)

szAx = g_z.add_subplot(312)
szAx.plot(time, s1z, 'ko-', linewidth=1.5, markevery=10)
szAx.plot(time, s2z, 'bo-', linewidth=1.5, markevery=10)
szAx.plot(time, s3z, 'ro-', linewidth=1.5, markevery=10)

szAx.set_ylim([0, 0.6])
szAx.set_ylabel(r'$\cos \theta = (\mathbf{s}_i \cdot \mathbf{z})$', fontsize=15)
szAx.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])
szAx.xaxis.set_minor_locator(AutoMinorLocator())
szAx.yaxis.set_minor_locator(AutoMinorLocator())
szAx.tick_params(which='major', length=6)
szAx.tick_params(which='minor', length=3)

wzAx = g_z.add_subplot(313)
wzAx.plot(time, wz, 'ko-', linewidth=1.5)

wzAx.set_xlabel('Time')
wzAx.set_ylabel(r'$\cos \theta = (\mathbf{\omega} \cdot \mathbf{z})$', fontsize=15)
wzAx.xaxis.set_minor_locator(AutoMinorLocator())
wzAx.yaxis.set_minor_locator(AutoMinorLocator())
wzAx.tick_params(which='major', length=6)
wzAx.tick_params(which='minor', length=3)
 

